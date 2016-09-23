from hydradna_config import *
import itertools
from hydradna.models import *
import warnings


class Reaction(object):
    @staticmethod
    def _anneal(template, primer, min_bases=MIN_BASES, threeprime=True):
        t = template.get_first()
        p = primer.get_first()
        if threeprime:
            t = t.find_last()
            p = p.find_last()
        visited = set()
        i = 0
        matches = []
        while t and t not in visited:
            next_method = lambda x: x.next()
            if threeprime:
                next_method = lambda x: x.prev()
            l = t._longest_match(p, next_method)
            if len(l) >= min_bases:
                matches.append((i, len(l)))
            visited.add(t)
            t = next_method(t)
            i += 1
        return matches

    @staticmethod
    def anneal_fiveprime(template, primer, min_bases=MIN_BASES):
        return Reaction._anneal(template, primer, min_bases=min_bases, threeprime=False)

    @staticmethod
    def anneal_threeprime(template, primer, min_bases=MIN_BASES):
        matches = Reaction._anneal(template, primer, min_bases=min_bases, threeprime=True)
        return [(len(template) - m[0], m[1]) for m in matches]

    @staticmethod
    def anneal_primer(template, primer, min_bases=MIN_BASES):
        rc_template = deepcopy(template).reverse_complement()
        fwd_matches = Reaction.anneal_threeprime(template, primer, min_bases=min_bases)
        # for f in fwd_matches:
        #     primer.cut
        rev_matches = Reaction.anneal_threeprime(rc_template, primer, min_bases=min_bases)

        def ca(matches):
            return [dict(pos=m[0], len=m[1], tm=Reaction.tm(primer.cut(len(primer) - m[1])[-1])) for m in matches]

        anneal = dict(F=ca(fwd_matches), R=ca(rev_matches))

        return anneal

    @staticmethod
    def polymerase(template, primer, min_bases, direction='forward'):
        products = []
        ann = Reaction.anneal_primer(template, primer, min_bases=min_bases)
        matches = None
        if direction=='forward':
            matches = ann['F']
        elif direction=='reverse':
            matches = ann['R']
        for match in matches:
                # New Product
                new_product = None

                # Overhang
                overhang = None
                try:
                    overhang, anneal = primer.cut(len(primer) - match['len'])
                except:
                    pass

                # Find binding position
                pos = match['pos'] - match['len']
                if direction=='forward':
                    # If forward...
                    new_product = template.cut(pos)[-1]
                    # Fuse overhang
                    if overhang is not None:
                        overhang.fuse(new_product)
                elif direction == 'reverse':
                    # If reverse...
                    new_product = deepcopy(template).reverse_complement().cut(pos)[-1]
                    # Fuse overhang
                    if overhang is not None:
                        overhang.fuse(new_product)
                    new_product.reverse_complement()
                products.append(new_product)
        return products, matches

    @staticmethod
    def _pcr(template, p1, p2, min_bases):
        def _polymerase_helper(template, primers, direction='forward'):
            x = []
            y = []
            for p in primers:
                products, matches = Reaction.polymerase(template, p, min_bases, direction=direction)
                x += products
                y += matches
            return zip(x, y)

        products = []
        matches = []
        for fwd_product, fwd_match in _polymerase_helper(template, [p1, p2], direction='forward'):
            for product, rev_match in _polymerase_helper(fwd_product, [p1, p2], direction='reverse'):
                products.append(product)
                matches.append((fwd_match, rev_match))
        return products, matches

    @staticmethod
    def pcr_analysis(template, p1, p2, min_bases):
        products, matches = Reaction._pcr(template, p1, p2, min_bases)
        report = ''


    @staticmethod
    def pcr(template, p1, p2, min_bases=MIN_BASES):
        products, matches = Reaction._pcr(template, p1, p2, min_bases)
        # rename products
        for p, m in zip(products, matches):
            print m
            p.name = '{} PCR[{}, {}]'.format(template.name, m[0]['pos'], len(template) - m[1]['pos'])
            p.description = """
                            PCR Product:\n
                            Input:\n
                            \tfwd_primer: {} {}\n""".format(str(p1), str(m[0]))
        return products

    @staticmethod
    def assembly_cycles(fragments, max_homology, min_homology):
        def complete_match(m):
            return m[0][0] == m[0][1]

        def less_than_max_homology(m):
            return m[0][0] < max_homology

        def one_match(m):
            return len(m) == 1

        def pass_conditions(m):
            return one_match(m) \
                   and less_than_max_homology(m) \
                   and complete_match(m)

        fragments += [deepcopy(f).reverse_complement() for f in fragments[1:]]
        pairs = itertools.product(fragments, fragments[:]) #itertools.permutations(fragments, 2)
        graph = []
        for pair in pairs:
            x = pair
            match = Reaction.anneal_threeprime(*pair, min_bases=min_homology)
            if match and pass_conditions(match):
                graph.append(pair)
        cycles = Utilities.Graph.find_cycles(graph)
        return cycles

    # TODO: make a linear assembly
    @staticmethod
    def cyclic_assembly(fragments, max_homology=MAX_GIBSON_HOMOLOGY, min_homology=MIN_BASES):
        cycles = Reaction.assembly_cycles(fragments, max_homology=max_homology, min_homology=min_homology)
        if len(cycles) > 1:
            warnings.warn("More than one assembly found.")
        products = []
        for cy in cycles:
            # Copy fragments in cycle
            cy = deepcopy(cy)

            # Pair fragments for cyclic assembly
            x1 = cy
            x2 = x1[1:] + x1[:1]
            pairs = zip(x1, x2)

            # Cut 5' ends of fragments according to homology
            for right, left in pairs:
                match = Reaction.anneal_threeprime(right, left)[0]
                nt = right.get()[match[0]]
                p = nt.cut_prev()
                right.first = nt

                # Fixing homology features
                homology1 = Sequence(first_nt=p)

                nt = left.get()[len(left) - match[0] - 1]
                n = nt.cut_next()
                left.first = nt
                homology2 = Sequence(first_nt=n)

                # add features for homology1
                if homology1.first is not None:
                    if not len(homology1) == len(homology2):
                        raise Exception("Homologies for assembly are different lengths.")

                for n1 in homology1.get():
                    for n2 in homology2.get():
                        n1.copy_features_from(n1)

                # fuse homology to left
                left.fuse(homology1)
            # Fuse all fragments
            for p in pairs:
                right, left = p
                left.fuse(right)

            products.append(x1[0])
        return products

    @staticmethod
    def end_homology(left, right):
        match = Reaction.anneal_threeprime(right, left)[0]
        nt = right.get()[match[0]]
        right.cut

    @staticmethod
    def tm(seq):
        if not (isinstance(seq, Sequence) or isinstance(seq, str)):
            raise TypeError("Cannot calculate tm of {} type for {}".format(str(type(seq)), seq))
        s = str(seq).lower()
        gc = s.count('g') + s.count('c')
        return 64.9 + 41 * (gc - 16.4) / len(s)