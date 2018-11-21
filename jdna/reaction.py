"""
Simulate molecular reactions
"""

import itertools
import warnings
from collections import defaultdict
from copy import copy
from jdna.sequence import Sequence
from jdna.graph import Graph

# Default values
MIN_BASES = 13
MAX_GIBSON_HOMOLOGY = 60


class Reaction(object):

    @staticmethod
    def combine_dnas(*parts):
        for p in parts:
            p.create_feature(p.name, 'misc', 0, len(p) - 1)
        last = None
        for p in parts:
            if last is None:
                last = p
                continue
            last = last + p
        return last


    #
    # # TODO: what does 'anneal' supposed to return???
    # @staticmethod
    # def _anneal(template, primer, min_bases=MIN_BASES, threeprime=True):
    #     t = template.head
    #     p = primer.head
    #     if threeprime:
    #         t = t.find_last()
    #         p = p.find_last()
    #     visited = set()
    #     i = 0
    #     matches = []
    #     while t and t not in visited:
    #         next_method = lambda x: next(x)
    #         if threeprime:
    #             next_method = lambda x: x.prev()
    #         l = t._longest_match(p, next_method)
    #         if len(l) >= min_bases:
    #             matches.append((i, len(l)))
    #         visited.add(t)
    #         t = next_method(t)
    #         i += 1
    #     return matches
    #
    # @staticmethod
    # def anneal_fiveprime(template, primer, min_bases=MIN_BASES):
    #     return Reaction._anneal(template, primer, min_bases=min_bases, threeprime=False)
    #
    # @staticmethod
    # def anneal_threeprime(template, primer, min_bases=MIN_BASES):
    #     matches = Reaction._anneal(template, primer, min_bases=min_bases, threeprime=True)
    #     return [(len(template) - m[0], m[1]) for m in matches]
    #
    # @staticmethod
    # def anneal_primer(template, primer, min_bases=MIN_BASES):
    #     fwd_matches = Reaction.anneal_threeprime(template, primer, min_bases=min_bases)
    #     # for f in fwd_matches:
    #     #     primer.cut
    #     template.reverse_complement()
    #     rev_matches = Reaction.anneal_threeprime(template, primer, min_bases=min_bases)
    #     template.reverse_complement()
    #
    #     def ca(matches):
    #         return [dict(primer=primer, pos=m[0], len=m[1], tm=Reaction.tm(primer.cut(len(primer) - m[1])[-1])) for m in
    #                 matches]
    #
    #     anneal = dict(F=ca(fwd_matches), R=ca(rev_matches))
    #
    #     return anneal

    @staticmethod
    def pcr(template, p1, p2, min_bases=MIN_BASES):
        ann1 = Reaction.anneal_primer(template, p1, min_bases=MIN_BASES)
        ann2 = Reaction.anneal_primer(template, p2, min_bases=MIN_BASES)

        f = ann1['F'] + ann2['F']
        r = ann1['R'] + ann2['R']

        pairs = itertools.product(f, r)
        products = []
        for p0, p1 in pairs:
            # make a new product
            product = copy(template)
            nts = product.nodes
            start_index = p0['pos'] - p0['len']
            end_index = len(template) - p1['pos'] + p1['len'] - 1
            s = nts[start_index]
            e = nts[end_index]
            s.cut_prev()
            e.cut_next()

            # if end is not in the product, continue
            product.head = s
            if e not in product.nodes:
                continue

            # get overhangs of primers
            fwd_primer = copy(p0['primer'])
            rev_primer = copy(p1['primer'])
            o1 = fwd_primer.get(len(fwd_primer) - p0['len']).cut_prev()
            o2 = rev_primer.get(len(rev_primer) - p1['len']).cut_prev()
            rev_primer.reverse_complement()

            # fuse overhangs
            if o1 is not None:
                Sequence(first=o1).fuse(product)
            if o2 is not None:
                product.fuse(Sequence(first=o2))
            products.append(product)
        return products

    @staticmethod
    def get_homology_graph(fragment_list, max_homology, min_homology):
        def complete_match(m):
            return m[0][0] == m[0][1]

        def less_than_max_homology(m):
            return m[0][0] <= max_homology

        def one_match(m):
            return len(m) == 1

        def pass_conditions(m):
            return one_match(m) \
                   and less_than_max_homology(m) \
                   and complete_match(m)

        fragments = fragment_list[:]

        # examine reversed sequences as well
        for f in fragments[1:]:
            reversed = copy(f).reverse_complement()
            if reversed.name is None:
                reversed.name = ''
            reversed.name = reversed.name + '(reversed)'
            fragments.append(reversed)

        fragment_to_id = {f: i for i, f in enumerate(fragments)}

        # examine all pairs
        pairs = itertools.product(fragments, fragments[:])  # itertools.permutations(fragments, 2)
        
        graph = defaultdict(list)
        match_graph = defaultdict(list)
        for p0, p1 in pairs:
            # TODO: add possibility of fragment to circularize on itself
            if p0 == p1:
                s = str(p0)
                if s[:min_homology] == s[-min_homology:]:
                    raise Exception("Fragment self circularized.")
            match = Reaction.anneal_threeprime(template=p1, primer=p0, min_bases=min_homology)
            # TODO: what is match? why is it a list??
            if match and pass_conditions(match) and match[0][0] < len(p1) and match[0][1] < len(p0):
                left, right = fragment_to_id[p0], fragment_to_id[p1]
                graph[left].append(right)
                match_graph[left].append(match)
        key = fragments
        return graph, key, match_graph

    # TODO: fix this method
    @staticmethod
    def homology_report(fragment_list, max_homology=MAX_GIBSON_HOMOLOGY, min_homology=MIN_BASES):
        fragment_list = [copy(f) for f in fragment_list]
        graph, fragments, match_graph = Reaction.get_homology_graph(fragment_list, max_homology, min_homology)
        cyclic_assemblies = Graph.find_cycles(graph)
        linear_assemblies = Graph.find_linear(graph)
        report = {
            'fragments': fragments,
            'interaction graph': graph,
            'homology graph': match_graph,
            'cyclic_paths': cyclic_assemblies,
            'linear_paths': linear_assemblies,
        }
        return report

    @staticmethod
    def print_homology_report(hr):
        c = hr['cyclic_paths']
        l = hr['linear_paths']
        f = hr['fragments']
        h = hr['homology graph']
        ig = hr['interaction graph']
        print(('Cyclic Assemblies: {}'.format(len(c))))
        for i, a in enumerate(c[::-1]):
            print(('\tAssemblyGraph {}'.format(i)))
            for n in a:
                print(('\t\t{} {}'.format(f[n].name, h[n])))
        print(('Linear Assemblies: {}'.format(len(l))))
        for i, a in enumerate(l[::-1]):
            print(('\tAssemblyGraph {}'.format(i)))
            for n in a:
                print(('\t\t{} {}'.format(f[n].name, h[n])))

    @staticmethod
    def cyclic_assembly(fragments, max_homology=MAX_GIBSON_HOMOLOGY, min_homology=MIN_BASES):
        return Reaction.homology_assembly(fragments, True, max_homology=max_homology, min_homology=min_homology)

    @staticmethod
    def linear_assembly(fragments, max_homology=MAX_GIBSON_HOMOLOGY, min_homology=MIN_BASES):
        return Reaction.homology_assembly(fragments, False, max_homology=max_homology, min_homology=min_homology)

    @staticmethod
    def homology_assembly(fragment_list, cyclic, max_homology=MAX_GIBSON_HOMOLOGY, min_homology=MIN_BASES):
        fragment_list = [copy(f) for f in fragment_list]
        h_report = Reaction.homology_report(fragment_list, max_homology=max_homology, min_homology=min_homology)
        fragments = h_report['fragments']
        paths = h_report['cyclic_paths']
        if not cyclic:
            paths = h_report['linear_paths']
        assembly = [[fragments[x] for x in y] for y in paths]
        if len(assembly) > 1:
            warnings.warn("More than one assembly found.")
        elif len(assembly) == 0:
            Reaction.print_homology_report(h_report)
            Exception("AssemblyGraph failed.")
        products = []
        for cyno, cy in enumerate(assembly):
            # Copy fragments in cycle
            cy = [copy(x) for x in cy]

            # Pair fragments for cyclic assembly
            x1 = cy[:-1]
            x2 = cy[1:]
            if cyclic:
                x1 = cy[:]
                x2 = cy[1:] + cy[:-1]
                # x2 = x1[:-1] + x1[1:]
            pairs = list(zip(x1, x2))
            # Cut 5' ends of fragments according to homology
            # overlap_info = []

            for left, right in pairs:
                # Find homology region and cut
                matches = Reaction.anneal_threeprime(right, left)
                if len(matches) == 0:
                    raise Exception(
                        "The fragment \"{}\" does not anneal to fragment \"{}\"".format(left.name, right.name))

                match = matches[0]

                right_nt = right.get(match[0])

                homology1_head_nt = right_nt.cut_prev()
                homology1 = Sequence(first=homology1_head_nt)
                right.head = right_nt

                left_nt = left.get(len(left) - match[0] - 1)
                homology2_head_nt = left_nt.cut_next()
                left.head = left_nt
                homology2 = Sequence(first=homology2_head_nt)

                # add features for homology1
                if homology1.head is not None:
                    if not len(homology1) == len(homology2):
                        raise Exception("Homologies for assembly are different lengths.")

                # copy features while removing overlapping ones
                for n1 in homology1.nodes:
                    for n2 in homology2.nodes:
                        n1.copy_features_from(n2)

                # overlap_info = (0, len(left), len(homology1), Reaction.tm(homology1), left.name)

                # fuse homology to left
                left.fuse(homology1)

            # Fuse all fragments
            for p in pairs:
                left, right = p
                left.fuse(right)

            # Annotate new product
            cyprod = x1[0]
            label = 'Cyclic: '
            if not cyclic:
                label = 'Linear'
            cyprod.name = '{} Product: {}'.format(label, cyno)

            products.append(cyprod)
        return products

    # def circularize_by_homology(self, seq, MAX_HOMOLOGY=MAX_GIBSON_HOMOLOGY, min_homology=MIN_BASES):
    #     seq_str = str(seq)
    #     seq_copy = copy(seq).chop_off_threeprime(MAX_HOMOLOGY)
    #     m = Reaction.anneal_threeprime(seq_copy, seq)
    #     if len(m) == 0:
    #         raise Exception("No homology found to circularize DNA.")
    #     elif len(m) > 1:
    #         raise Exception("More than one homology site found to circularize DNA.")
    #     match = m[0]
    #     if not match[0] == match[1]:
    #         raise Exception("Homology site is not a complete match. Cannot circularize DNA.")
    #     seq.chop_off_threeprime(len(seq) - match[0] - 1)
    #     seq.make_cyclic()
    #     if not str(seq) == seq_str:
    #         raise Exception("Error occured while attempting circularization by homology. Sequences do not match.")
    #     return seq

    # TODO: finish overlap extension pcr
    @staticmethod
    def overlap_extension_pcr(fragment_list, primer1, primer2, max_homology=MAX_GIBSON_HOMOLOGY,
                              min_homology=MIN_BASES):
        fragment_list = [copy(f) for f in fragment_list]
        graph, fragments, match_graph = Reaction.get_homology_graph(fragment_list, max_homology, min_homology)
        print((type(graph)))
        print((list(graph.keys())))
        print(graph)
        linear_assemblies = Graph.find_linear(graph)
        print(linear_assemblies)

    # TODO: add test modules for this
    @staticmethod
    def end_homology(left, right, copy=True):
        if copy:
            left = copy(left)
            right = copy(right)

        # Find homology region and cut
        match = Reaction.anneal_threeprime(right, left)[0]
        nt = right.get(match[0])
        p = nt.cut_prev()

        # modify right fragment
        right.head = nt

        # define homology 1
        homology1 = Sequence(first=p)

        # define homology 2
        nt = left.get(len(left) - match[0] - 1)
        n = nt.cut_next()
        left.head = nt
        homology2 = Sequence(first=n)

        # add features for homology1
        if homology1.head is not None:
            if not len(homology1) == len(homology2):
                raise Exception("Homologies for assembly are different lengths.")

        # copy features while removing overlapping ones
        for n1 in homology1.nodes:
            for n2 in homology2.nodes:
                n1.copy_features_from(n2)

        # overlap_info = (0, len(left), len(homology1), Reaction.tm(homology1), left.name)

        # fuse homology to left
        left.fuse(homology1)
        left.fuse(right)

        left.get_first()
        return left, homology2

    @staticmethod
    def tm(seq):
        if not (isinstance(seq, Sequence) or isinstance(seq, str)):
            raise TypeError("Cannot calculate tm of {} type for {}".format(str(type(seq)), seq))
        s = str(seq).lower()
        gc = s.count('g') + s.count('c')
        return 64.9 + 41 * (gc - 16.4) / len(s)
