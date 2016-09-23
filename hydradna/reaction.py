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
    def polymerase(template, primer, direction='forward'):
        products = []
        ann = Reaction.anneal_primer(template, primer)
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
        return products

    @staticmethod
    def pcr(template, p1, p2):
        def _polymerase_helper(template, primers, direction='forward'):
            x = []
            for p in primers:
                x += Reaction.polymerase(template, p, direction=direction)
            return x

        products = []
        fwd_products = _polymerase_helper(template, [p1, p2], direction='forward')
        for fwd_product in fwd_products:
            rev_products = _polymerase_helper(fwd_product, [p1, p2], direction='reverse')
            products += rev_products
        return products

    # @staticmethod
    # def pcr_analysis(template, p1, p2, min_bases=MIN_BASES):
    #     products, fwd_report, rev_report = Reaction._pcr(template, p1, p2, min_bases=min_bases)
    #     primer_dict = {p1: 1, p2: 2}
    #     print products
    #     print "***  PCR  ***"
    #     print "-" * 10
    #     print "INPUT:"
    #     print ">Template length: {}".format(len(template))
    #     print ">Primer 1: {}".format(p1)
    #     print ">Primer 2: {}".format(p2)
    #     print "-"*10
    #     print "OUTPUT:"
    #     print "Num Products: {}".format(len(products))
    #     for i, z in enumerate(zip(products, fwd_report, rev_report)):
    #         p, f, r = z
    #         print
    #         print "\tProduct {}".format(i)
    #         print "\t Length: {}".format(len(p))
    #         print "\t\tTm_anneal\tTm_whole\tAnnealLen\t5'Pos"
    #         print "\tF_Bind (p{}):".format(primer_dict[f['primer']]), \
    #             '\t\t'.join([str(f[x]) for x in \
    #                          ['tm_anneal', 'tm_whole', 'anneal_length', 'anneal_pos']])
    #         print "\tR_Bind (p{}):".format(primer_dict[r['primer']]), \
    #             '\t\t'.join([str(r[x]) for x in \
    #                          ['tm_anneal', 'tm_whole', 'anneal_length', 'anneal_pos']])


    # @staticmethod
    # def pcr(template, p1, p2, min_bases=MIN_BASES):
    #
    #     # Annealing locations for p1
    #     ann1 = Reaction.anneal_primer(template, p1, min_bases=min_bases)
    #
    #     # Annealing locations for p2
    #     ann2 = Reaction.anneal_primer(template, p2, min_bases=min_bases)
    #
    #     # Group fwd and rev matches
    #     fwd_matches = ann1['F'] + ann2['F']
    #     rev_matches = ann1['R'] + ann2['R']
    #
    #     # All possible paired matches
    #     pairs = itertools.product(fwd_matches, rev_matches)
    #     print list(pairs)
    #     pairs = itertools.product(fwd_matches, rev_matches)
    #     products = []
    #     for fwd, rev in pairs:
    #         # Copy the template
    #         template_copy = deepcopy(template)
    #
    #         # Find the start binding position
    #         s = fwd['pos'] - fwd['len']
    #         start = template_copy.get()[s]
    #
    #         # Cut the start
    #         start.cut_prev()
    #
    #         # Simulate the polymerase
    #         fwd_product =
    #
    #         # Find binding positions on the template
    #         s = fwd['pos'] - fwd['len']
    #         e = len(template_copy) - rev['pos'] + rev['len'] - 1
    #         start = template_copy.get()[s]
    #         end = template_copy.get()[e]
    #
    #         # Move forward from start to end
    #         # Simulate the polymerase
    #         postart.fwd(stop_link=end)
    #
    #         start.cut_prev()
    #         end.cut_next()
    #         product = Sequence(first_nt=start)
    #
    #         # Fuse the overhangs
    #         overhang1 = p1.cut(len(p1) - fwd['len'])[0]
    #         overhang2 = p2.cut(len(p2) - rev['len'])[0]
    #         product.insert(overhang1, 0)
    #         product.insert(overhang2.reverse_complement(), len(product))
    #         products.append(product)
    #     return products

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

            # Pair fragments for assembly
            x1 = cy
            x2 = x1[1:] + x1[:1]
            pairs = zip(x1, x2)

            # Cut 5' ends of fragments according to homology
            for right, left in pairs:
                match = Reaction.anneal_threeprime(right, left)[0]
                nt = right.get()[match[0]]
                nt.cut_prev()
                right.first = nt

            # Fuse all fragments
            for p in pairs:
                right, left = p
                left.fuse(right)

            products.append(x1[0])
        return products

    @staticmethod
    def tm(seq):
        if not (isinstance(seq, Sequence) or isinstance(seq, str)):
            raise TypeError("Cannot calculate tm of {} type for {}".format(str(type(seq)), seq))
        s = str(seq).lower()
        gc = s.count('g') + s.count('c')
        return 64.9 + 41 * (gc - 16.4) / len(s)