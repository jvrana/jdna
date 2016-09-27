from collections import defaultdict
import itertools
import warnings
from copy import copy, deepcopy
import random
import re

# Default values
MIN_BASES = 13
MAX_GIBSON_HOMOLOGY = 60


class Link(object):
    """
    A link represents a single piece of data in a linked sequence.
    Each link can be connected to another link through the _next
    and _prev pointers (representing backbone bonds).

    Available methods:
       get_next
       get_prev
       add_next
       add_prev
       cut_next
       cut_prev
       remove_next
       remove_prev
       find_first
       find_last
       fwd
       rev
    """

    def __init__(self, data):
        self.data = data
        self.__next = None
        self.__prev = None

    def next(self):
        return self.__next

    def prev(self):
        return self.__prev

    def add_next(self, data):
        new_link = Link(data)
        self.set_next(new_link)
        return new_link

    def add_prev(self, data):
        new_link = Link(data)
        self.set_prev(new_link)
        return new_link

    def cut_next(self):
        next_link = self.next()
        if next_link is not None:
            next_link.__assign_prev(None)
        self.__assign_next(None)
        return next_link

    def cut_prev(self):
        prev_link = self.prev()
        if prev_link is not None:
            prev_link.__assign_next(None)
        self.__assign_prev(None)
        return prev_link

    def break_connections(self):
        self.set_next(None)
        self.set_prev(None)

    def remove(self):
        next_link = self.next()
        prev_link = self.prev()
        if next_link is not None:
            next_link.set_prev(prev_link)
        if prev_link is not None:
            prev_link.set_next(next_link)
        self.break_connections()
        return

    def swap(self):
        temp = self.__next
        self.__next = self.__prev
        self.__prev = temp

    def __assign_next(self, link):
        self.__next = link

    def __assign_prev(self, link):
        self.__prev = link

    def set_next(self, link):
        if link is not None:
            link.__assign_prev(self)
        self.__assign_next(link)

    def set_prev(self, link):
        if link is not None:
            link.__assign_next(self)
        self.__assign_prev(link)

    def make_cyclic(self):
        if not self.is_cyclic():
            first = self.find_first()
            last = self.find_last()
            last.set_next(first)

    def is_cyclic(self):
        visited = set()
        curr = self
        while curr:
            if curr in visited:
                return True
            visited.add(curr)
            curr = curr.next()
        return False

    def _propogate(self, next_method, stop=None, stop_criteria=None):
        visited = []
        nxt = self
        while nxt is not None:
            if stop_criteria is not None and stop_criteria(nxt):
                break
            visited.append(nxt)
            if nxt is stop:
                break
            nxt = next_method(nxt)
            if nxt is visited[0]:
                break
        return visited

    def fwd(self, stop_link=None, stop_criteria=None):
        return self._propogate(
            lambda x: x.next(),
            stop=stop_link,
            stop_criteria=stop_criteria)

    def rev(self, stop_link=None, stop_criteria=None):
        return self._propogate(
            lambda x: x.prev(),
            stop=stop_link,
            stop_criteria=stop_criteria)

    def find_first(self):
        links = self.rev()
        return links[-1]

    def find_last(self):
        links = self.fwd()
        return links[-1]

    def _longest_match(self, y, next_method):
        x1 = self
        x2 = y
        longest_match = []
        while x1 and x2:
            if x1.equivalent(x2):
                longest_match.append((x1, x2))
                x1 = next_method(x1)
                x2 = next_method(x2)
            else:
                break
        return longest_match

    def _complete_match(self, y, next_method):
        l = self._longest_match(y, next_method)
        if not l:
            return False
        t1, t2 = self._longest_match(y, next_method)[-1]
        return not (next_method(t1) and next_method(t2))

    def complete_match_fwd(self, y):
        return self._complete_match(y, lambda x: x.next())

    def complete_match_rev(self, y):
        return self._complete_match(y, lambda x: x.prev())

    def equivalent(self, other):
        return self.data == other.data

    def __copy__(self):
        copied = type(self)(self.data)
        return copied

    def __deepcopy__(self, memo):
        raise NotImplementedError("copy.deepcopy not implemented with class"\
                "{}. Use copy.copy instead.".format(self.__class__.__name__))

    def __repr__(self):
        return str(self)

    def __str__(self):
        return str(self.data)


class DoubleLinkedList(object):

    def __init__(self, first=None, sequence=None):
        if sequence is not None:
            self.initialize(sequence)
        elif first is not None:
            self.first = first

    def initialize(self, sequence):
        self.first = Link(sequence[0])
        current = self.first
        for d in sequence[1:]:
            new = Link(d)
            current.set_next(new)
            current = new

    def get_first(self):
        if self.is_cyclic():
            return self.first
        first = self.first.find_first()
        self.first = first
        return self.first

    def set_first(self, link):
        self.first = link

    def make_cyclic(self):
        return self.get_first().make_cyclic()

    def is_cyclic(self):
        visited = set()
        curr = self.first
        while curr:
            if curr in visited:
                return True
            visited.add(curr)
            curr = curr.next()
        return False

    def linearize(self, i=0):
        this_i = self.get()[i]
        this_i.cut_prev()
        return this_i

    def get(self):
        return self.get_first().fwd()

    def cut(self, i, cut_prev=True):
        if isinstance(i, tuple):
            i = list(i)
        if isinstance(i, int):
            i = [i]
        # Special case in which i == len
        i = list(set(i))
        if len(self) in i and cut_prev:
            i.remove(len(self))
            if self.is_cyclic():
                i.append(0)
        i = list(set(i))
        i.sort()
        self._inbounds(i)
        self_copy = copy(self)
        all_links = self_copy.get()
        for cut_loc in i:
            link = all_links[cut_loc]
            if cut_prev:
                link.cut_prev()
            else:
                link.cut_next()
        return DoubleLinkedList._group_links(all_links)

    #TODO: Speed up with set
    @staticmethod
    def _group_links(links):
        unique_first_links = []
        for link in links:
            first_link = link.find_first()
            if first_link not in unique_first_links:
                unique_first_links.append(first_link)
        return [DoubleLinkedList(first=l) for l in unique_first_links]

    def insert(self, linkedlist, i, copy_insertion=True):
        if i == len(self.get()):
            pass
        else:
            self._inbounds(i)
        if linkedlist.is_cyclic():
            raise TypeError("Cannot insert a cyclic sequence")
        if copy_insertion:
            linkedlist = copy(linkedlist)
        #TODO: This copies the insertion sequence, you want that?
        if i == len(self.get()):
            loc2 = None
            loc1 = self.get()[i-1]
        else:
            loc2 = self.get()[i]
            loc1 = loc2.prev()
        first = linkedlist.get()[0]
        last = linkedlist.get()[-1]
        first.set_prev(loc1)
        last.set_next(loc2)
        if i == 0:  # Special case in which user inserts sequence in front of their sequence; they probably intend to re-index it
            self.first = first
        return self

    def remove(self, i):
        self._inbounds(i)
        to_be_removed = self.get()[i]
        new_first = self.get_first()
        if i == 0:
            new_first = new_first.next()
        to_be_removed.remove()
        self.first = new_first
        return

    def reindex(self, i):
        self._inbounds(i)
        if not self.is_cyclic():
            raise TypeError("Cannot re-index a linear linked set")
        self.first = self.get()[i]

    def _inbounds(self, num):
        if isinstance(num, int):
            num = [num]
        for n in num:
            mn = 0
            mx = len(self.get()) - 1
            if n < 0 or n > mx:
                raise IndexError("Index {} out of acceptable bounds ({}, {})".format(n, mn, mx))

    def search_all(self, query):
        curr_link = self.get_first()
        q_link = query.get_first()
        i = 0
        found = []
        visited = set()
        while curr_link and curr_link not in visited:
            visited.add(curr_link)
            if curr_link.complete_match_fwd(q_link):
                found.append((i, curr_link))
            curr_link = curr_link.next()
            i += 1
        return found

    def reverse(self):
        for s in self.get():
            s.swap()
        if self.is_cyclic():
            self.reindex(1)
        return self

    def slice(self, i, j, fwd=True):
        links = self.get()
        start = links[i]
        stop = links[j]
        method = start.fwd
        if not fwd:
            method = start.rev
        sec_links = method(stop_link=stop)
        if sec_links[-1] is stop:
            return sec_links
        else:
            raise IndexError("Improper indices for linkedlist.")


    def __copy__(self):
        copied = type(self)(first=Link(''))
        copied.__dict__.update(self.__dict__)
        copied.first = copy(self.get_first())
        curr = copied.first
        for link in self.get()[1:]:
            copied_link = copy(link)
            curr.set_next(copied_link)
            curr = copied_link
        if self.is_cyclic():
            copied.make_cyclic()
        return copied

    def __deepcopy__(self, memo):
        raise NotImplementedError("copy.deepcopy not implemented with class" \
                              "{}. Use copy.copy instead.".format(self.__class__.__name__))

    def __reversed__(self):
        for s in self.get():
            s.swap()
        return self

    def __len__(self):
        return len(self.get())

    def __iter__(self):
        current = self.first
        while current is not None:
            yield current
            current = current.next()

    def __repr__(self):
        return str(self)

    def __str__(self):
        return ''.join(str(x) for x in self.get())


class Feature(object):

    def __init__(self, name, type='misc feature', strand=1, color=None):
        self.name = name
        self.type = type
        self.strand = 1
        self.length = None
        if color is None:
            color = random_color()
        self.color = color

    def __str__(self):
        return '{} {}'.format(self.name, self.type)

    def __repr__(self):
        return str(self)


def rgb_to_hex(r, g, b):
    def clamp(x):
        return max(0, min(x, 255))

    return "#{0:02x}{1:02x}{2:02x}".format(clamp(r), clamp(g), clamp(b))

def random_color():
    rgb = [int(random.random()*255) for x in range(3)]
    return rgb_to_hex(*rgb)


class Nucleotide(Link):

    base_pairing = dict(zip(
        ['a', 't', 'c', 'g', 'A', 'T', 'C', 'G'],
        ['t', 'a', 'g', 'c', 'T', 'A', 'G', 'C']
    ))

    def __init__(self, base):
        super(Nucleotide, self).__init__(base)
        self.features = {}

    def base(self):
        return self.data

    def equivalent(self, other):
        return str(self.base()).upper() == str(other.base()).upper()

    def to_complement(self):
        self.data = self.base_pairing[self.data]

    def set_next(self, nucleotide):
        super(Nucleotide, self).set_next(nucleotide)
        Nucleotide.fuse_features(self, nucleotide)

    def set_prev(self, nucleotide):
        super(Nucleotide, self).set_prev(nucleotide)
        Nucleotide.fuse_features(nucleotide, self)

    def cut_prev(self):
        return self._cut(cut_prev=True)

    def cut_next(self):
        return self._cut(cut_prev=False)

    def _cut(self, cut_prev=True):
        for f in self.features:
            self.split_features(split_prev=cut_prev)
        nxt = None
        if cut_prev:
            nxt = super(Nucleotide, self).cut_prev()
        else:
            nxt = super(Nucleotide, self).cut_next()
        return nxt

    def feature_fwd(self, feature):
        stop = lambda x: feature not in x.features
        return self._propogate(lambda x: x.next(), stop_criteria=stop)

    def feature_rev(self, feature):
        stop = lambda x: feature not in x.features
        return self._propogate(lambda x: x.prev(), stop_criteria=stop)

    def add_feature(self, feature, pos):
        self.features[feature] = pos
        return feature

    def replace_feature(self, old_feature, new_feature):
        self.features[new_feature] = self.features[old_feature]
        self.remove_feature(old_feature)

    def copy_features_from(self, other):
        for f in other.features:
            i = other.features[f]
            if f not in self.features:
                self.add_feature(f, i)
        self._remove_overlapping_features()

    def remove_feature(self, feature):
        del self.features[feature]

    def get_feature_span(self, feature):
        start = self.feature_rev(feature)[-1]
        end = self.feature_fwd(feature)[-1]
        return (start.features[feature], end.features[feature])

    # def update_feature_span(self, feature, delta_i):
    #     start = self.feature_rev(feature)[-1]
    #     for n in start.feature_fwd(feature):
    #         n.features[feature] += delta_i

    def _remove_overlapping_features(self):
        # type: () -> Nucleotide
        feature_pairs = itertools.combinations(self.features.keys(), 2)
        tobedel = set()
        for f1, f2 in feature_pairs:
            if f1.name == f2.name:
                tobedel.add(f2)
        for tob in tobedel:
            self.remove_feature(tob)

    @staticmethod
    def fuse_features(n1, n2):
        if n1 is None:
            return
        if n2 is None:
            return

        delset = set()

        for f1 in n1.features:
            for f2 in n2.features:
                f1_pos = n1.features[f1]
                f2_pos = n2.features[f2]
                f1_copy = copy(f1)
                # same name & consecutive position
                if f1 is f2:
                    continue
                if f1.name == f2.name and f1_pos + 1 == f2_pos:
                    delset.add((f1, f2, f1_copy))
        for f1, f2, f1_copy in delset:
            for n in n1.feature_rev(f1):
                try:
                    n.replace_feature(f1, f1_copy)
                except KeyError:
                    pass
            for n in n2.feature_fwd(f2):
                try:
                    n.replace_feature(f2, f1_copy)
                except KeyError:
                    pass


    def split_features(self, split_prev=True):
        x1 = self.prev()
        x2 = self
        if not split_prev:
            # then split_next
            x1 = self
            x2 = self.next()
        # If at the end, no splitting is necessary
        if x1 is None or x2 is None:
            return
        for f in x1.features:
            # If this feature spans
            if f in x2.features:
                # Grab the sequences for the split feature
                frag1 = x1.feature_rev(f)
                frag2 = x2.feature_fwd(f)

                # Make two copies of the feature
                f1 = copy(f)
                f2 = copy(f)

                # Swap original feature for copy
                for n in frag1:
                    n.replace_feature(f, f1)
                for n in frag2:
                    n.replace_feature(f, f2)

    def __copy__(self):
        copied = super(Nucleotide, self).__copy__()
        copied.features = deepcopy(self.features)
        return copied


class Sequence(DoubleLinkedList):

    def __init__(self, first=None, sequence=None, name='unknown'):
        super(Sequence, self).__init__(first=first, sequence=sequence)
        self.name = name
        self.description = ''

    def initialize(self, sequence):
        self.first = Nucleotide(sequence[0])
        current = self.first
        for d in sequence[1:]:
            new = Nucleotide(d)
            current.set_next(new)
            current = new

    def add_feature(self, i, j, feature, start_index=0):
        s = self.get()
        start = s[i]
        end = s[j]
        feature_nts = start.fwd(stop_link=end)
        if not feature_nts[-1] == end:
            raise IndexError("Feature index error")
        for i, n in enumerate(feature_nts):
            n.add_feature(feature, i+start_index)
        feature.length = len(feature_nts)
        return feature

    def _features_to_i(self):
        features = defaultdict(list)
        for i, x in enumerate(self.get()):
            for f in x.features:
                features[f].append((x, i))
        return features

    def get_features(self):
        features_to_nts = self._features_to_i()
        feature_info = {}
        for feature in features_to_nts:
            nt_to_i_list = features_to_nts[feature]
            nts, indices = zip(*nt_to_i_list)
            first = nts[0].feature_rev(feature)[-1]
            last = nts[0].feature_fwd(feature)[-1]
            nt_to_i = dict(zip(self.get(), range(len(self))))
            feature_range = (first.features[feature], last.features[feature])
            pos_ranges = (nt_to_i[first], nt_to_i[last])
            feature_info[feature] = [pos_ranges, feature_range]
        return feature_info

    def print_features(self):
        features = self.get_feature_pos()
        print "Features:",
        for f in features:
            e = f.end
            if e is None:
                e = ''
            span = '[{}...{}]'.format(f.start, e)
            if f.end == None and f.start == 0:
                span = ''
            print '{}{} {}'.format(f.name, span, f.type),
            ranges = list(Utilities.group_ranges(features[f]))
            feature_ranges = []
            for r in ranges:
                print '{}...{},'.format(r[0], r[-1]),
        print

    def find_feature(self, name=None):
        found = []
        for feature in self.get_features():
            if feature.name == name:
                found.append(feature)
        return found

    def create_feature(self, name, feature_type, start, end):
        f = Feature(name, feature_type)
        self.add_feature(start, end, f)
        return f

    def complement(self):
        curr = self.get_first()
        visited = set()
        while curr and curr not in visited:
            visited.add(curr)
            curr.to_complement()
            curr = curr.next()
        return self

    def reverse_complement(self):
        self.reverse()
        self.complement()
        return self

    def cut(self, i, cut_prev=True):
        fragments = super(Sequence, self).cut(i, cut_prev)
        fragments = [Sequence(first=f.get_first()) for f in fragments]
        return fragments

    def fuse(self, seq):
        f = self.get_first().find_last()
        l = seq.get_first()
        f.set_next(l)
        return self

    def __copy__(self):
        return super(Sequence, self).__copy__()

    def __add__(self, other):
        return copy(self).fuse(copy(other))

    def __repr__(self):
        return str(self)


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
        rc_template = copy(template).reverse_complement()
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
                    new_product = copy(template).reverse_complement().cut(pos)[-1]
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
    def pcr(template, p1, p2, min_bases=MIN_BASES):
        products, matches = Reaction._pcr(template, p1, p2, min_bases)
        # rename products
        for p, m in zip(products, matches):
            p.name = '{} PCR[{}, {}]'.format(template.name, m[0]['pos'], len(template) - m[1]['pos'])
            report = """*** PCR ***
                        INPUTS:
                        Template: {template}
                        P1: {primer1} {match1}
                        P2: {primer2} {match2}

                        OUTPUS:
                        Products: {numproducts}""".format(**dict(
                template=template.name,
                primer1=str(p1),
                match1=str(m[0]),
                primer2=str(p2),
                match2=str(m[1]),
                numproducts=len(products)
            ))
            p.description = report
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

        fragments += [copy(f).reverse_complement() for f in fragments[1:]]
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
            cy = copy(cy)

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
                homology1 = Sequence(first=p)

                nt = left.get()[len(left) - match[0] - 1]
                n = nt.cut_next()
                left.first = nt
                homology2 = Sequence(first=n)

                # add features for homology1
                if homology1.first is not None:
                    if not len(homology1) == len(homology2):
                        raise Exception("Homologies for assembly are different lengths.")

                # copy features while removing overlapping ones
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


class Utilities:


    class Graph:

        @staticmethod
        def find_cycles(graph, min_path_length=1):
            def find_new_cycles(path):
                start_node = path[0]
                next_node = None
                sub = []

                # visit each edge and each node of each edge
                for edge in graph:
                    node1, node2 = edge
                    if start_node == node1:
                        next_node = node2
                    if not visited(next_node, path):
                        # neighbor node not on path yet
                        sub = [next_node]
                        sub.extend(path)
                        # explore extended path
                        find_new_cycles(sub)
                    elif len(path) >= min_path_length and next_node == path[-1]:
                        # cycle found
                        p = invert(path)
                        p = rotate_to_smallest(p)
                        if is_new(p): # and isNew(inv):
                            cycles.append(p)

            def invert(pth):
                return rotate_to_smallest(pth[::-1])

            #  rotate cycle path such that it begins with the smallest node
            def rotate_to_smallest(path):
                n = path.index(min(path))
                return path[n:] + path[:n]

            def is_new(path):
                return path not in cycles

            def visited(node, path):
                return node in path

            cycles = []
            for e in graph:
                for n in e:
                    find_new_cycles([n])
            return cycles

    @staticmethod
    def group_ranges(lst):
        pos = (j - i for i, j in enumerate(lst))
        t = 0
        for i, els in itertools.groupby(pos):
            l = len(list(els))
            el = lst[t]
            t += l
            yield range(el, el + l)


class Convert:

    @staticmethod
    def from_benchling(data):

        def clean_data(dic):
            for key in dic:
                if isinstance(dic[key], basestring):
                    dic[key] = str(dic[key])

        clean_data(data)
        # Sequence, name
        seq = Sequence(sequence=data['bases'], name=data['name'])

        # Topology
        if data['circular']:
            seq.make_cyclic()

        def fix_end(end):
            end = end - 1
            if end < 0:
                end = len(seq) - 1
            return end

        for a in data['annotations']:
            # Feature name, type, strand, color

            name, span_start, span_end, span_length = Convert.parse_feature_name(a['name'])
            newf = Feature(name, a['type'], strand=a['strand'], color=a['color'])
            start_index = 0
            if span_start is not None:
                start_index = int(span_start)
            seq.add_feature(a['start'], fix_end(a['end']), newf, start_index=start_index)
            if span_length is not None:
                newf.length = int(span_length)
        return seq

    @staticmethod
    def to_benchling_json(seq):
        data = {}
        data['name'] = seq.name
        data['circular'] = seq.is_cyclic()
        data['annotations'] = []
        data['bases'] = str(seq)
        data['description'] = seq.description

        def add_annotation(f):
            data['annotations'].append(f)

        def fix_end(end):
            end += 1
            if end == len(seq):
                end = 0
            return end

        features = seq.get_features()
        for f in features:
            range, span = features[f]
            start, end = range
            span = '[{},{},{}]'.format(span[0], span[1], f.length)
            add_annotation(dict(
                name=str(f.name) + ' ' + str(span),
                type=f.type,
                start=start,
                end=fix_end(end),
                strand=f.strand,
                color=f.color
            ))
        return data

    @staticmethod
    def parse_feature_name(fname):
        g = re.search('(.+)(\[(\d+),(\d*),(\d*)\])|(.+)', fname)
        name, span, span_start, span_end, span_length, fullname = g.groups()
        if name is None:
            name = fullname
        return name, span_start, span_end, span_length