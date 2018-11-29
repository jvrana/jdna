"""
Simulate molecular reactions
"""

import itertools
from collections import namedtuple

import networkx as nx

from jdna.sequence import Sequence, SequenceFlags
from jdna.viewer import SequenceViewer
import primer3
from decimal import Decimal

class ReactionException(Exception):
    """Generic reaction exception"""


class PCRException(ReactionException):
    """Exception with pcr"""


BindingEvent = namedtuple("BindingEvent", ["template", "primer", "position"])


class Assembly(object):

    def __init__(self, templates, overhangs, cyclic):
        self.templates = [t.copy() for t in templates]
        self.overhangs = [o.copy() for o in overhangs]
        self.cyclic = cyclic

    def tms(self):
        return [round(o.tm(), 2) for o in self.overhangs]

    @property
    def product(self):
        product = Sequence()
        if not self.cyclic:
            zipped = zip(self.templates, self.overhangs)
        else:
            zipped = zip(self.overhangs, self.templates)
        for s1, s2 in zipped:
            product += s1 + s2
        if self.cyclic:
            product.circularize()
        return product

    def view(self):
        seqs = []
        positions = [0]
        pos = 0
        overhangs = self.overhangs[:]
        for o in overhangs:
            o.annotate(0, len(o)-1, "Tm: {}°C".format(round(o.tm(), 1)))
        if self.cyclic:
            overhangs.append(overhangs[0])
        else:
            overhangs.append(Sequence())
        for i, t in enumerate(self.templates):
            seq = overhangs[0] + t + overhangs[i + 1]
            seqs.append(overhangs[i] + t + overhangs[i + 1])
            pos += len(seq) - len(overhangs[i + 1])
            positions.append(pos)

        aligned_seqs = []
        for p, s in zip(positions, seqs):
            aligned_seqs.append(Sequence('-' * p) + s)

        mx = max([len(s) for s in aligned_seqs])
        for i, s in enumerate(aligned_seqs):
            diff = mx - len(s)
            aligned_seqs[i] = aligned_seqs[i] + Sequence('-' * diff)

        viewer = SequenceViewer(aligned_seqs, sequence_labels=['({i}) {index}'.format(i=i, index="{index}") for i in
                                                               range(len(aligned_seqs))])
        for seq in aligned_seqs:
            Sequence._apply_features_to_view(seq, viewer)
        viewer.metadata['Num Fragments'] = len(self.templates)
        viewer.metadata['Overhang Tms (°C)'] = ', '.join([str(x) for x in self.tms()])
        viewer.metadata['Overhang Lengths (bp)'] = ', '.join([str(len(x)) for x in self.overhangs])
        viewer.metadata['Overhang ΔG'] = self.format_float_array(self.deltaGs())
        viewer.metadata['Overhang ΔG (hairpin)'] = self.format_array(self.deltaG_hairpins())
        viewer.metadata['Competing ΔG'] = self.format_float_array(self.competing_deltaGs())
        viewer.metadata['Length'] = "{}bp".format(len(self.product))
        return viewer

    @staticmethod
    def format_array(arr):
        return ', '.join([str(x) for x in arr])

    @staticmethod
    def format_float_array(arr):
        return ', '.join([
            "{:.2e}".format(Decimal(float(x))) for x in arr
        ])

    def deltaG_hairpins(self):
        gs = []
        for o in self.overhangs:
            fwd = primer3.calcHairpin(str(o).upper()).dg
            rev = primer3.calcHairpin(str(o.copy().rc()).upper()).dg
            gs.append((fwd, rev))
        return gs

    def deltaGs(self):
        return [primer3.calcHeterodimer(str(o).upper(), str(o.copy().reverse_complement()).upper()).dg for o in self.overhangs]

    def competing_deltaGs(self):
        gs = []
        for o in self.overhangs:
            total = 0
            total += primer3.calcHairpin(str(o).upper()).dg
            total += primer3.calcHairpin(str(o.copy().rc()).upper()).dg
            others = self.overhangs[:]
            others.remove(o)
            for other in others:
                dg1 = primer3.calcHeterodimer(str(o).upper(), str(other.copy().rc()).upper()).dg
                dg2 = primer3.calcHeterodimer(str(o).upper(), str(other).upper()).dg
                total += dg1
                total += dg2
            gs.append(total)
        return gs

    def print(self):
        self.view().print()

    def __str__(self):
        return str(self.view())


class Reaction(object):
    MIN_BASES = 13

    @classmethod
    def pcr(cls, template, p1, p2, min_bases=MIN_BASES):
        bindings = list(template.anneal(p1, min_bases=min_bases))
        bindings += template.anneal(p2, min_bases=min_bases)

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
                "Some primers did not bind. Number of forward bindings: {}. Number of rev bindings: {}".format(
                    len(forward_bindings), len(reverse_bindings)
                ))
        products = []
        for fwd, rev in itertools.product(forward_bindings, reverse_bindings):
            overhang1 = fwd.five_prime_overhang
            overhang2 = rev.five_prime_overhang
            print(fwd.span)
            print(rev.span)
            print(template.index_of(fwd.start))
            print(template.index_of(rev.end))
            amplified = template.new_slice(fwd.start, rev.end)
            print(overhang1)
            print(overhang2)
            products.append(overhang1 + amplified + overhang2)
        return products

    @classmethod
    def anneal_sequences(cls, sequences, min_bases=Sequence.DEFAULTS.MIN_ANNEAL_BASES):
        # pairs = itertools.product(sequences, sequences + [s.copy().reverse_complement() for s in sequences[1:]])
        pairs = itertools.product(sequences, repeat=2)
        bindings = []
        for s1, s2 in pairs:
            for binding in s1.dsanneal(s2, min_bases=min_bases):
                if not (binding.span[0] == 0 and binding.span[-1] == len(s1) - 1):
                    bindings.append(BindingEvent(s1, s2, binding))
        return bindings

    @classmethod
    def interaction_graph(cls, sequences, min_bases=Sequence.DEFAULTS.MIN_ANNEAL_BASES, bind_reverse_complement=False,
                          depth=None):
        G = nx.DiGraph()
        if bind_reverse_complement:
            sequences = [s.copy() for s in sequences] + [s.copy().reverse_complement() for s in sequences[:]]
        for s in sequences:
            G.add_node(s.global_id, sequence=s)
        seq_pairs = itertools.product(sequences, repeat=2)
        for s1, s2 in seq_pairs:
            for binding in s1.anneal_forward(s2, min_bases=min_bases, depth=depth):
                if not (binding.span[0] == 0 and binding.span[-1] == len(s1) - 1):
                    G.add_edge(s2.global_id, s1.global_id, template=s1, primer=s2, binding=binding)
        return G

    @classmethod
    def linear_paths(cls, G):
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
        paths = []
        for path in nx.simple_cycles(G.copy()):
            paths.append(path)
        return paths

    @classmethod
    def path_to_edge_data(cls, G, path, cyclic):
        if cyclic:
            path_pairs = zip(path, path[1:] + [path[0]])
        else:
            path_pairs = zip(path[:-1], path[1:])
        edges = []
        for n1, n2 in path_pairs:
            edges.append(G.edges[n1, n2])
        return edges

    @classmethod
    def path_to_assembly(self, path, G, cyclic=False):
        prev_template_end = None

        node_pairs = []
        overhangs = []

        for i, data in enumerate(self.path_to_edge_data(G, path, cyclic)):
            binding = data['binding']

            # use the previous template and this query should refer to the same sequence
            node_pairs.append((prev_template_end, binding.query_start.prev()))

            # the next segment start is the template start
            prev_template_end = binding.end.next()
            overhangs.append(binding.template_anneal)

        if not cyclic:
            node_pairs.append((prev_template_end, None))
            overhangs.append(Sequence())
        else:
            pass
            node_pairs[0] = (prev_template_end, node_pairs[0][1])
            overhangs = [overhangs[-1]] + overhangs[:-1]
        amplified_sequences = [Sequence.new_slice(*pair) for pair in node_pairs]

        return Assembly(amplified_sequences, overhangs, cyclic)

    @classmethod
    def linear_assemblies(cls, sequences, min_bases=Sequence.DEFAULTS.MIN_ANNEAL_BASES, depth=None):
        G = cls.interaction_graph(sequences, bind_reverse_complement=True, min_bases=min_bases, depth=depth)
        linear_paths = cls.linear_paths(G)
        if not linear_paths:
            return []
        assemblies = []
        for path in linear_paths:
            assembly = cls.path_to_assembly(path, G)
            assemblies.append(assembly)
        return assemblies

    @classmethod
    def cyclic_assemblies(cls, sequences, min_bases=Sequence.DEFAULTS.MIN_ANNEAL_BASES, depth=None):
        G = cls.interaction_graph(sequences, bind_reverse_complement=True, min_bases=min_bases, depth=depth)

        nodes = list(G.nodes)
        fwd = nodes[:len(sequences)]
        rev = nodes[len(sequences):]
        fpriority = list(range(len(fwd)))
        rpriority = list(range(len(fwd), len(fwd) * 2))
        priority_rank = fpriority + rpriority[1:] + [rpriority[0]]
        node_priority = dict(zip(fwd + rev, priority_rank))

        cyclic_paths = cls.cyclic_paths(G)
        if not cyclic_paths:
            return []
        assemblies = []
        cyclic_paths = sorted(cyclic_paths, key=lambda x: min([node_priority[n] for n in x]))
        for path in cyclic_paths:
            ranks = [node_priority[n] for n in path]
            mn = ranks.index(min(ranks))
            cyclic_path = path[mn:] + path[:mn]
            assembly = cls.path_to_assembly(cyclic_path, G, cyclic=True)
            assemblies.append(assembly)
        return assemblies
