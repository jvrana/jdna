"""Simulate molecular reactions."""
import itertools
from collections import namedtuple
from decimal import Decimal

import networkx as nx
import primer3

from jdna.sequence import Sequence
from jdna.sequence import SequenceFlags
from jdna.viewer import SequenceViewer


class ReactionException(Exception):
    """Generic reaction exception."""


class PCRException(ReactionException):
    """Exception with pcr."""


BindingEvent = namedtuple("BindingEvent", ["template", "primer", "position"])


class Assembly:
    """Stores a DNA assembly.

    Assembly can be printed to display relevant assembly information:

    .. code::

        assembly.print()

    ::

        > "Unnamed" (320bp)
          Cyclic: True
          Num Fragments: 4
          Overhang Tms (°C): 53.65, 53.79, 54.96, 56.37
          Junction Lengths (bp): 20, 20, 20, 20
          Junction ΔG: -1.91e+4, -1.93e+4, -1.98e+4, -2.05e+4
          Junction ΔG (hairpin): (0.0, 232.50050218543765), (203.8205043739763,
          -84.3344956260189), (0.0, -60.91449781456504),
          (-1026.7134934374844, -622.0224934343823)
          Competing ΔG: -2.18e+4, -2.11e+4, -1.92e+4, -2.01e+4
          Product length (bp): 300

                  -----Tm: 53.7°C-----
        (0) 0     TTGACGACAGTGGCTATCCCCTGTTGCTAGGCACGGTGATATATCAGCCC
        (1) 0     --------------------------------------------------
        (2) 0     --------------------------------------------------
        (3) 0     --------------------------------------------------

                                           -----Tm: 53.8°C-----
        (0) 50    ATAGGGCCGGATCAACGTAGTTTGATCGTAGCTGTTCCGTCAGTC-----
        (1) 50    -------------------------TCGTAGCTGTTCCGTCAGTCAGTGC
        (2) 50    --------------------------------------------------
        (3) 50    --------------------------------------------------

        (0) 100   --------------------------------------------------
        (1) 100   CGATGCCAGTCTACTGCTTTTCGCCCAGGGGGACACCTGACACTATGTTA
        (2) 100   --------------------------------------------------
        (3) 100   --------------------------------------------------

                  -----Tm: 55.0°C-----
        (0) 150   --------------------------------------------------
        (1) 150   TGCCAAGCTCACGTTTCACA------------------------------
        (2) 150   TGCCAAGCTCACGTTTCACAAGGCAGATTAAATAAGACTGAAACATTTGG
        (3) 150   --------------------------------------------------

                                           -----Tm: 56.4°C-----
        (0) 200   --------------------------------------------------
        (1) 200   --------------------------------------------------
        (2) 200   TGGAGGCGCGAGTCTTGCCCTTAGGCAGACCTGTTCGAGGGCGAA-----
        (3) 200   -------------------------CAGACCTGTTCGAGGGCGAAAAACG

        (0) 250   --------------------------------------------------
        (1) 250   --------------------------------------------------
        (2) 250   --------------------------------------------------
        (3) 250   CCGAAGTCACATAGCGGCATGTTAGGTAGCTGTACACCCGTACTGCTACA
    """

    def __init__(self, templates, junctions, cyclic):
        """Assembly constructor.

        :param templates: list of templates
        :type templates: list
        :param junctions: list of junctions
        :type junctions: list
        :param cyclic: whether assembly is cyclic
        :type cyclic: bool
        """
        self._templates = tuple([t.copy() for t in templates])
        self._junctions = tuple([o.copy() for o in junctions])
        self._cyclic = cyclic

    @property
    def templates(self):
        return tuple(t.copy() for t in self._templates)

    @property
    def junctions(self):
        return tuple(j.copy() for j in self._junctions)

    @property
    def cyclic(self):
        return self._cyclic

    def tms(self):
        """Return the tms of the junctions."""
        return [round(o.tm(), 2) for o in self.junctions]

    @property
    def product(self):
        product = Sequence()
        product.name = "Unnamed assembly"
        if not self.cyclic:
            zipped = zip(self.templates, self.junctions)
        else:
            zipped = zip(self.junctions, self.templates)
        for s1, s2 in zipped:
            product += s1 + s2
        if self.cyclic:
            product.circularize()
        return product

    def view(self, *args, **kwargs):
        seqs = []
        positions = [0]
        pos = 0
        junctions = list(self.junctions)
        for jxn in junctions:
            jxn.annotate(None, None, "Tm: {}°C".format(round(jxn.tm(), 1)))
        if self.cyclic:
            junctions.append(junctions[0])
        else:
            junctions = [Sequence()] + list(junctions)
        for i, t in enumerate(self.templates):
            seq = junctions[i] + t + junctions[i + 1]
            seqs.append(seq)
            pos += len(seq) - len(junctions[i + 1])
            positions.append(pos)

        aligned_seqs = []
        for p, s in zip(positions, seqs):
            aligned_seqs.append(Sequence("-" * p) + s)

        mx = max([len(s) for s in aligned_seqs])
        for i, s in enumerate(aligned_seqs):
            diff = mx - len(s)
            aligned_seqs[i] = aligned_seqs[i] + Sequence("-" * diff)

        labels = [
            "({i}) {index}".format(i=i, index="{index}")
            for i in range(len(aligned_seqs))
        ]
        viewer = SequenceViewer(
            aligned_seqs,
            *args,
            sequence_labels=labels,
            foreground_colors="RANDOM",
            **kwargs
        )
        for seq in aligned_seqs:
            Sequence._apply_features_to_view(seq, viewer)
        viewer.metadata["Cyclic"] = self.cyclic
        viewer.metadata["Num Fragments"] = len(self.templates)
        viewer.metadata["Fragment Names"] = "\n\t" + "\n\t".join(
            ["{} - {}".format(i, t.name) for i, t in enumerate(self.templates)]
        )
        viewer.metadata["Overhang Tms (°C)"] = ", ".join([str(x) for x in self.tms()])
        viewer.metadata["Junction Lengths (bp)"] = ", ".join(
            [str(len(x)) for x in self.junctions]
        )
        viewer.metadata["Junction ΔG"] = self.format_float_array(self.deltaGs())
        viewer.metadata["Junction ΔG (hairpin)"] = self.format_array(
            self.deltaG_hairpins()
        )
        viewer.metadata["Competing ΔG"] = self.format_float_array(
            self.competing_deltaGs()
        )
        viewer.metadata["Product length (bp)"] = "{}".format(len(self.product))
        return viewer

    @staticmethod
    def format_array(arr):
        return ", ".join([str(x) for x in arr])

    @staticmethod
    def format_float_array(arr):
        return ", ".join(["{:.2e}".format(Decimal(float(x))) for x in arr])

    def deltaG_hairpins(self):
        gs = []
        for jxn in self.junctions:
            fwd = self._hairpin(jxn).dg
            rev = self._hairpin(jxn.copy().rc()).dg
            gs.append((fwd, rev))
        return gs

    def _heterodimer(self, s1, s2):
        return primer3.calcHeterodimer(str(s1[:60]).upper(), str(s2[:60]).upper())

    def _hairpin(self, s):
        return primer3.calcHairpin(str(s[:60]).upper())

    def deltaGs(self):
        overhangs = [o[:60] for o in self.junctions]
        return [self._heterodimer(o, o.copy().rc()).dg for o in overhangs]

    def competing_deltaGs(self):
        gs = []
        jxns = list(self.junctions)
        for jxn in jxns:
            total = 0
            total += self._hairpin(jxn).dg
            total += self._hairpin(jxn.copy().rc()).dg
            others = jxns[:]
            others.remove(jxn)
            for other in others:
                dg1 = self._heterodimer(jxn, other.copy().rc()).dg
                dg2 = self._heterodimer(jxn, other.copy()).dg
                total += dg1
                total += dg2
            gs.append(total)
        return gs

    def print(self, *args, **kwargs):
        """
        Print the assembly.

        ::

            > "Unnamed" (320bp)
                  Num Fragments: 4
                  Overhang Tms (°C): 61.49, 58.37, 54.72, 45.91
                  Junction Lengths (bp): 20, 20, 20, 20
                  Junction ΔG: -2.30e+4, -2.18e+4, -1.97e+4, -1.54e+4
                  Junction ΔG (hairpin): (-1542.5069956260231, -1787.4319956260224),
                  (-1832.7404978145642, -708.4919956260237), (0.0, 0.0), (0.0, 0.0)
                  Competing ΔG: -3.39e+4, -2.43e+4, -2.19e+4, -1.33e+4
                  Length: 300bp

                                                                                                     |<Tm: 58.4°C
                          -----Tm: 61.5°C-----                                                       ----------
                (0) 0     TGCTCGGCGAAGGGCCTGACAAGACCACTTCGTACCCTTTGTAACGTGACTTCTGTGAATCGATGGCGGATGTTGCTTGTGACGC
                (1) 0     ---------------------------------------------------------------------------CTTGTGACGC
                (2) 0     -------------------------------------------------------------------------------------
                (3) 0     -------------------------------------------------------------------------------------

                          |<Tm: 58.4°C
                          ----------                                                       -----Tm: 54.7°C-----
                (0) 85    ACCCGTAGCG---------------------------------------------------------------------------
                (1) 85    ACCCGTAGCGACATAACCGCCTTATTCCCACTTGTTTGGGGAGGCAAGTTTTGTAGTAGCCATTCTAATCCCGTTTTCTCCGCCG
                (2) 85    -----------------------------------------------------------------TAATCCCGTTTTCTCCGCCG
                (3) 85    -------------------------------------------------------------------------------------

                                                                                 -----Tm: 45.9°C-----
                (0) 170   -------------------------------------------------------------------------------------
                (1) 170   -------------------------------------------------------------------------------------
                (2) 170   CTTGTTCAGCTTGATGTGTTGGAACAGCAAGTTTATTTTAGGGTGTCAGGAGGGGGAGGTCCTCTAGTATTAAGC----------
                (3) 170   -------------------------------------------------------GAGGTCCTCTAGTATTAAGCACCGGCCAAG

                                                                       -----Tm: 61.5°C-----
                (0) 255   -----------------------------------------------------------------
                (1) 255   -----------------------------------------------------------------
                (2) 255   -----------------------------------------------------------------
                (3) 255   GCGGCACGGACGCCACTAGAAGGGAGGCGTTCATAGCGAGTCCTTTGCTCGGCGAAGGGCCTGAC

        """
        self.view(*args, **kwargs).print()

    def __str__(self):
        return str(self.view())


class Reaction:
    MIN_BASES = 13

    @classmethod
    def pcr(cls, template, primers, min_bases=MIN_BASES):
        """Make a pcr product from a template and primers.

        :param template: template
        :type template: Sequence
        :param primers: primers
        :type primers: list | tuple | iterable
        :param min_bases: minimum number of bases for primer annealing
        :type min_bases: int
        :return: list of Sequence instances
        :rtype: list
        """
        bindings = []
        for p in primers:
            bindings += template.anneal(p, min_bases=min_bases)

        forward_bindings = []
        reverse_bindings = []
        for bind in bindings:
            if bind.direction == SequenceFlags.FORWARD:
                forward_bindings.append(bind)
            elif bind.direction == SequenceFlags.REVERSE:
                reverse_bindings.append(bind)
            else:
                raise Exception("Direction {} not recognized".format(bind.direction))

        if not forward_bindings or not reverse_bindings:
            raise PCRException(
                "Some primers did not bind. Number of forward bindings: {}. Number of "
                "rev bindings: {}".format(len(forward_bindings), len(reverse_bindings))
            )
        products = []
        for fwd, rev in itertools.product(forward_bindings, reverse_bindings):
            overhang1 = fwd.five_prime_overhang
            overhang2 = rev.five_prime_overhang
            amplified = template.new_slice(fwd.start, rev.end)
            products.append(
                overhang1 + amplified + overhang2.copy().reverse_complement()
            )
        return products

    @classmethod
    def interaction_graph(
        cls,
        sequences,
        min_bases=Sequence.DEFAULTS.MIN_ANNEAL_BASES,
        bind_reverse_complement=False,
        max_bases=None,
        only_ends=False,
    ):
        """Make an interaction graph from a list of sequences.

        :param sequences: list of sequences
        :type sequences: list
        :param min_bases: minimum number of bases to bind
        :type min_bases: int
        :param bind_reverse_complement: whether to reverse complement each sequence and
                check interactions
        :type bind_reverse_complement: int
        :param max_bases: max number of bases to search
        :type max_bases: int
        :param only_ends: whether to only consider end-to-end interactions
        :type only_ends: bool
        :return:
        :rtype:
        """
        G = nx.DiGraph()
        sequences = [s.copy() for s in sequences]
        if bind_reverse_complement:
            reverse_complement = []
            for s in sequences:
                new = s.copy().reverse_complement()
                new.name = s.name + "(-1)"
                reverse_complement.append(new)
            sequences += reverse_complement
        for s in sequences:
            G.add_node(s.global_id, sequence=s)
        seq_pairs = itertools.product(sequences, repeat=2)
        for s1, s2 in seq_pairs:
            for binding in s1.anneal_forward(s2, min_bases=min_bases, depth=max_bases):
                if only_ends:
                    if binding.span[0] != 0 or binding.query_span[-1] != len(s2) - 1:
                        break
                if not (binding.span[0] == 0 and binding.span[-1] == len(s1) - 1):
                    G.add_edge(
                        s2.global_id,
                        s1.global_id,
                        template=s1,
                        primer=s2,
                        binding=binding,
                    )
        return G

    @classmethod
    def interaction_report(cls, G):
        rows = []
        for n1, n2 in G.edges:
            edge = G.edges[n1, n2]
            binding = edge["binding"]
            template = edge["template"]
            primer = edge["primer"]
            template_info = "{name} (id={gid} ({start},{end}/{length})".format(
                name=template.name,
                gid=template.global_id,
                start=binding.span[0],
                end=binding.span[1],
                length=len(template) - 1,
            )
            primer_info = "{name} (id={gid} ({start},{end}/{length})".format(
                name=primer.name,
                gid=primer.global_id,
                start=binding.query_span[0],
                end=binding.query_span[1],
                length=len(primer) - 1,
            )
            rows.append("{} -> {}".format(template_info, primer_info))
        print("format: primer -> template")
        print("\n".join(rows))

    @classmethod
    def linear_paths(cls, G):
        """Find linear paths from an interaction graph."""
        paths = []
        subgraph = G.subgraph(G.nodes).copy()
        for _ in nx.simple_cycles(G):
            return []
        while subgraph.nodes:
            path = nx.dag_longest_path(subgraph)
            if len(path) < 2:
                break
            subgraph.remove_nodes_from(path)
            paths.append(path)
        return paths

    @classmethod
    def cyclic_paths(cls, G):
        """Find cyclic paths from an interaction graph."""
        paths = []
        for path in nx.simple_cycles(G.copy()):
            paths.append(path)
        return paths

    @classmethod
    def _path_to_edge_data(cls, G, path, cyclic):
        if cyclic:
            path_pairs = zip(path, path[1:] + [path[0]])
        else:
            path_pairs = zip(path[:-1], path[1:])
        edges = []
        for n1, n2 in path_pairs:
            edges.append(G.edges[n1, n2])
        return edges

    @classmethod
    def _path_to_assembly(self, path, G, cyclic=False):
        prev_template_end = None

        node_pairs = []
        overhangs = []

        for i, data in enumerate(self._path_to_edge_data(G, path, cyclic)):
            binding = data["binding"]

            # use the previous template and this query should refer to the same sequence
            node_pairs.append(
                (prev_template_end, binding.query_start.prev(), data["primer"].name)
            )

            # the next segment start is the template start
            prev_template_end = binding.end.next()
            overhangs.append(binding.template_anneal)

        if not cyclic:
            node_pairs.append((prev_template_end, None, data["template"].name))
            overhangs.append(Sequence())
        else:
            pass
            node_pairs[0] = (prev_template_end, node_pairs[0][1], node_pairs[0][2])
            overhangs = [overhangs[-1]] + overhangs[:-1]

        amplified_sequences = []
        for n1, n2, name in node_pairs:
            sequence = Sequence.new_slice(n1, n2)
            sequence.name = name
            amplified_sequences.append(sequence)

        return Assembly(amplified_sequences, overhangs, cyclic)

    @classmethod
    def linear_assemblies(
        cls, sequences, min_bases=Sequence.DEFAULTS.MIN_ANNEAL_BASES, max_bases=None
    ):
        """Finds all unique, longest linear assemblies. If cyclic assemblies
        are found, returns empty list.

        :param sequences: list of sequences
        :type sequences: list
        :param min_bases: minimum junction length
        :type min_bases: int
        :param max_bases: maximum junction length. Smaller values will result in faster
                    search, but may miss valid
                            assemblies.
        :return: list of Assembly instances. Sequence can be accessed using
                    `assembly.product`
        :rtype: list
        """
        G = cls.interaction_graph(
            sequences,
            bind_reverse_complement=True,
            min_bases=min_bases,
            max_bases=max_bases,
            only_ends=True,
        )
        linear_paths = cls.linear_paths(G)
        if not linear_paths:
            return []
        assemblies = []
        for path in linear_paths:
            assembly = cls._path_to_assembly(path, G)
            assemblies.append(assembly)
        return assemblies

    @classmethod
    def cyclic_assemblies(
        cls, sequences, min_bases=Sequence.DEFAULTS.MIN_ANNEAL_BASES, max_bases=None
    ):
        """Finds all unique, longest cyclic assemblies. If no cyclic assemblies
        found, returns empty list.

        :param sequences: list of sequences
        :type sequences: list
        :param min_bases: minimum junction length
        :type min_bases: int
        :param max_bases: maximum junction length. Smaller values will result in faster
                search, but may miss valid
                            assemblies.
        :type max_bases: int
        :return: list of Assembly instances. Sequence can be accessed using
                    `assembly.product`
        :rtype: list
        """
        G = cls.interaction_graph(
            sequences,
            bind_reverse_complement=True,
            min_bases=min_bases,
            max_bases=max_bases,
            only_ends=True,
        )

        nodes = list(G.nodes)
        fwd = nodes[: len(sequences)]
        rev = nodes[len(sequences) :]
        fpriority = list(range(len(fwd)))
        rpriority = list(range(len(fwd), len(fwd) * 2))
        priority_rank = fpriority + rpriority[1:] + [rpriority[0]]
        node_priority = dict(zip(fwd + rev, priority_rank))

        cyclic_paths = cls.cyclic_paths(G)
        if not cyclic_paths:
            return []
        assemblies = []
        cyclic_paths = sorted(
            cyclic_paths, key=lambda x: min([node_priority[n] for n in x])
        )
        for path in cyclic_paths:
            ranks = [node_priority[n] for n in path]
            mn = ranks.index(min(ranks))
            cyclic_path = path[mn:] + path[:mn]
            assembly = cls._path_to_assembly(cyclic_path, G, cyclic=True)
            assemblies.append(assembly)
        return assemblies
