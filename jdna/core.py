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

    def __next__(self):
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
        next_link = next(self)
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
        next_link = next(self)
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
            curr = next(curr)
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
            lambda x: next(x),
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
        return self._complete_match(y, lambda x: next(x))

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
            curr = next(curr)
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
        cut_links = []
        for cut_loc in i:
            link = all_links[cut_loc]
            c = None
            if cut_prev:
                c = link.cut_prev()
                if c is not None:
                    cut_links.append(c)
                cut_links.append(link)
            else:
                cut_links.append(link)
                c = link.cut_next()
                if c is not None:
                    cut_links.append(c)
        return DoubleLinkedList._group_links(cut_links)

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
            new_first = next(new_first)
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
            curr_link = next(curr_link)
            i += 1
        return found

    def reverse(self):
        for s in self.get():
            s.swap()
        if self.is_cyclic():
            self.reindex(1)
        return self

    # def slice(self, i, j, fwd=True):
    #     links = self.get()
    #     start = links[i]
    #     stop = links[j]
    #     method = start.fwd
    #     if not fwd:
    #         method = start.rev
    #     sec_links = method(stop_link=stop)
    #     if sec_links[-1] is stop:
    #         return sec_links
    #     else:
    #         raise IndexError("Improper indices for linkedlist.")


    def __copy__(self):
        copied = type(self)(sequence='X')
        copied.__dict__.update(self.__dict__)
        copied.initialize(str(self))
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
            current = next(current)

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

    base_pairing = dict(list(zip(
        ['a', 't', 'c', 'g', 'A', 'T', 'C', 'G'],
        ['t', 'a', 'g', 'c', 'T', 'A', 'G', 'C']
    )))

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
        self.cut_next()
        super(Nucleotide, self).set_next(nucleotide)
        Nucleotide.fuse_features(self, nucleotide)

    def set_prev(self, nucleotide):
        self.cut_prev()
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
        return self._propogate(lambda x: next(x), stop_criteria=stop)

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
        feature_pairs = itertools.combinations(list(self.features.keys()), 2)
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
            x2 = next(self)
        # If at the end, no splitting is necessary
        if x1 is None or x2 is None:
            return
        for f in x1.features:
            # If this feature spans
            if f in x2.features:
                # Grab the sequences for the split feature
                frag1 = x1.feature_rev(f)
                frag2 = x2.feature_fwd(f)

                # check if its a cyclic feature
                if x2 in frag1:
                    continue
                if x1 in frag2:
                    continue

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
            nts, indices = list(zip(*nt_to_i_list))
            first = nts[0].feature_rev(feature)[-1]
            last = nts[0].feature_fwd(feature)[-1]
            if next(last) == nts[0]:
                first = nts[0]
            nt_to_i = dict(list(zip(self.get(), list(range(len(self))))))
            feature_range = (first.features[feature], last.features[feature])
            pos_ranges = (nt_to_i[first], nt_to_i[last])
            feature_info[feature] = [pos_ranges, feature_range]
        return feature_info

    def print_features(self):
        raise NotImplementedError()
        # features = self.get_feature_pos()
        # print "Features:",
        # for f in features:
        #     e = f.end
        #     if e is None:
        #         e = ''
        #     span = '[{}...{}]'.format(f.start, e)
        #     if f.end == None and f.start == 0:
        #         span = ''
        #     print '{}{} {}'.format(f.name, span, f.type),
        #     ranges = list(Utilities.group_ranges(features[f]))
        #     feature_ranges = []
        #     for r in ranges:
        #         print '{}...{},'.format(r[0], r[-1]),
        # print

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
            curr = next(curr)
        return self

    def reverse_complement(self):
        self.reverse()
        self.complement()
        return self

    def cut(self, i, cut_prev=True):
        fragments = super(Sequence, self).cut(i, cut_prev)
        fragments = [Sequence(first=f.get_first()) for f in fragments]
        return fragments

    def chop_off_fiveprime(self, i):
        if self.is_cyclic():
            raise IndexError('Cannot chop a cyclic sequence.')
        return self.cut(i)[-1]

    def chop_off_threeprime(self, i):
        if self.is_cyclic():
            raise IndexError('Cannot chop a cyclic sequence.')
        return self.cut(i, cut_prev=False)[0]

    def fuse(self, seq):
        f = self.get_first().find_last()
        l = seq.get_first()
        f.set_next(l)
        return self

    def __copy__(self):
        copied = super(Sequence, self).__copy__()
        features = self.get_features()
        for f in features:
            pos, span = features[f]
            copied.add_feature(pos[0], pos[1], copy(f), start_index=span[0])
        return copied

    def __add__(self, other):
        return copy(self).fuse(copy(other))

    def __repr__(self):
        return str(self)


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
            next_method = lambda x: next(x)
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
        fwd_matches = Reaction.anneal_threeprime(template, primer, min_bases=min_bases)
        # for f in fwd_matches:
        #     primer.cut
        template.reverse_complement()
        rev_matches = Reaction.anneal_threeprime(template, primer, min_bases=min_bases)
        template.reverse_complement()
        def ca(matches):
            return [dict(primer=primer, pos=m[0], len=m[1], tm=Reaction.tm(primer.cut(len(primer) - m[1])[-1])) for m in matches]

        anneal = dict(F=ca(fwd_matches), R=ca(rev_matches))

        return anneal

    # @staticmethod
    # def polymerase(template, primer, min_bases, direction='forward'):
    #     products = []
    #     ann = Reaction.anneal_primer(template, primer, min_bases=min_bases)
    #     matches = None
    #     if direction=='forward':
    #         matches = ann['F']
    #     elif direction=='reverse':
    #         matches = ann['R']
    #     for match in matches:
    #             # New Product
    #             new_product = None
    #
    #             # Overhang
    #             overhang = None
    #             try:
    #                 overhang, anneal = primer.cut(len(primer) - match['len'])
    #             except:
    #                 pass
    #
    #             # Find binding position
    #             pos = match['pos'] - match['len']
    #             if direction=='forward':
    #                 # If forward...
    #                 new_product = template.cut(pos)[-1]
    #                 # Fuse overhang
    #                 if overhang is not None:
    #                     overhang.fuse(new_product)
    #             elif direction == 'reverse':
    #                 # If reverse...
    #                 new_product = copy(template).reverse_complement().cut(pos)[-1]
    #                 # Fuse overhang
    #                 if overhang is not None:
    #                     overhang.fuse(new_product)
    #                 new_product.reverse_complement()
    #             products.append(new_product)
    #     return products, matches

    # @staticmethod
    # def _pcr(template, p1, p2, min_bases):
    #     def _polymerase_helper(template, primers, direction='forward'):
    #         x = []
    #         y = []
    #         for p in primers:
    #             products, matches = Reaction.polymerase(template, p, min_bases, direction=direction)
    #             x += products
    #             y += matches
    #         return zip(x, y)
    #
    #     products = []
    #     matches = []
    #     for fwd_product, fwd_match in _polymerase_helper(template, [p1, p2], direction='forward'):
    #         for product, rev_match in _polymerase_helper(fwd_product, [p1, p2], direction='reverse'):
    #             products.append(product)
    #             matches.append((fwd_match, rev_match))
    #     return products, matches

    @staticmethod
    def pcr(template, p1, p2, min_bases=MIN_BASES):
        ann1 = Reaction.anneal_primer(template, p1, min_bases=MIN_BASES)
        ann2 = Reaction.anneal_primer(template, p2, min_bases=MIN_BASES)

        f = ann1['F'] + ann2['F']
        r = ann1['R'] + ann2['R']

        pairs = itertools.product(f, r)
        products = []
        for pair in pairs:
            # make a new product
            product = copy(template)
            nts = product.get()
            start_index = pair[0]['pos'] - pair[0]['len']
            end_index = len(template) - pair[1]['pos'] + pair[1]['len'] - 1
            s = nts[start_index]
            e = nts[end_index]
            s.cut_prev()
            e.cut_next()

            # if end is not in the product, continue
            product.first = s
            if e not in product.get():
                continue

            # get overhangs of primers
            fwd_primer = copy(pair[0]['primer'])
            rev_primer = copy(pair[1]['primer'])
            o1 = fwd_primer.get()[len(fwd_primer) - pair[0]['len']].cut_prev()
            o2 = rev_primer.get()[len(rev_primer) - pair[1]['len']].cut_prev()
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
        for f in fragments[1:]:
            reversed = copy(f).reverse_complement()
            reversed.name = reversed.name + '(reversed)'
            fragments.append(reversed)
        fragment_to_id = {}
        for i, f in enumerate(fragments):
            fragment_to_id[f] = i
        pairs = itertools.product(fragments, fragments[:]) #itertools.permutations(fragments, 2)
        graph = defaultdict(list)
        match_graph = defaultdict(list)
        for pair in pairs:
            # TODO: add possibility of fragment to circularize on itself
            if pair[0] == pair[1]:
                s = str(pair[0])
                if s[:min_homology] == s[-min_homology:]:
                    raise Exception("Fragment self circularized.")
            match = Reaction.anneal_threeprime(*pair, min_bases=min_homology)
            if match and pass_conditions(match):
                left, right = fragment_to_id[pair[0]], fragment_to_id[pair[1]]
                graph[left].append(right)
                match_graph[left].append(match)
        key = fragments
        return graph, key, match_graph

    # TODO: fix this method
    @staticmethod
    def homology_report(fragment_list, max_homology=MAX_GIBSON_HOMOLOGY, min_homology=MIN_BASES):
        fragment_list = [copy(f) for f in fragment_list]
        graph, fragments, match_graph = Reaction.get_homology_graph(fragment_list, max_homology, min_homology)
        cyclic_assemblies = Utilities.Graph.find_cycles(graph)
        linear_assemblies = Utilities.Graph.find_linear(graph)
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
            x1 = cy
            x2 = x1[1:] + x1[:1]
            pairs = list(zip(x1, x2))
            # Cut 5' ends of fragments according to homology
            # overlap_info = []

            for right, left in pairs:
                # Find homology region and cut
                match = Reaction.anneal_threeprime(right, left)[0]

                nt = right.get()[match[0]]
                p = nt.cut_prev()

                # modify right fragment
                right.first = nt

                # define homology sequences
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
                        n1.copy_features_from(n2)

                # overlap_info = (0, len(left), len(homology1), Reaction.tm(homology1), left.name)

                # fuse homology to left
                left.fuse(homology1)

            # Fuse all fragments
            for p in pairs:
                right, left = p
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
    def overlap_extension_pcr(fragment_list, primer1, primer2, max_homology=MAX_GIBSON_HOMOLOGY, min_homology=MIN_BASES):
        fragment_list = [copy(f) for f in fragment_list]
        graph, fragments, match_graph = Reaction.get_homology_graph(fragment_list, max_homology, min_homology)
        print((type(graph)))
        print((list(graph.keys())))
        print(graph)
        linear_assemblies = Utilities.Graph.find_linear(graph)
        print(linear_assemblies)

    # TODO: add test modules for this
    @staticmethod
    def end_homology(left, right, copy=True):
        if copy:
            left = copy(left)
            right = copy(right)

        # Find homology region and cut
        match = Reaction.anneal_threeprime(right, left)[0]
        nt = right.get()[match[0]]
        p = nt.cut_prev()

        # modify right fragment
        right.first = nt

        # define homology 1
        homology1 = Sequence(first=p)

        # define homology 2
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


class Utilities:


    class Graph:

        @staticmethod
        def find_cycles(graph, min_path_length=2):

            unique_cycles = []
            for g in graph:
                cycles = Utilities.Graph.find_cyclic_paths(graph, [g], min_path_length=min_path_length)
                for cy in cycles:
                    cy = Utilities.rotate_to_smallest(cy)
                    if cy not in unique_cycles:
                        unique_cycles.append(cy)
            return unique_cycles

        @staticmethod
        def find_linear(graph, min_path_length=2):
            unique_paths = []
            for g in graph:
                c = True
                for p in unique_paths:
                    if g in p:
                        c = False
                        break
                if c:
                    unique_paths += Utilities.Graph.find_linear_paths(graph, [g], min_path_length=min_path_length)
            return unique_paths

        @staticmethod
        def find_linear_paths(graph, path, min_path_length=2):
            if path[-1] not in graph and len(path) >= min_path_length:
                return [path]
            paths = []
            for node in graph[path[-1]]:
                if node not in path:
                    newpaths = Utilities.Graph.find_linear_paths(graph, path[:] + [node], min_path_length=min_path_length)
                    paths += newpaths
            return paths

        @staticmethod
        def find_cyclic_paths(graph, path, min_path_length=2):
            if len(path) >= min_path_length and path[-1] == path[0]:
                return [path[:-1]]
            paths = []
            if path[-1] not in graph:
                return []
            for node in graph[path[-1]]:
                if node not in path[1:]:
                    newpaths = Utilities.Graph.find_cyclic_paths(graph, path[:] + [node], min_path_length=min_path_length)
                    paths += newpaths
            return paths

    @staticmethod
    def rotate_to_smallest(path):
        n = path.index(min(path))
        return path[n:] + path[:n]

    @staticmethod
    def group_ranges(lst):
        pos = (j - i for i, j in enumerate(lst))
        t = 0
        for i, els in itertools.groupby(pos):
            l = len(list(els))
            el = lst[t]
            t += l
            yield list(range(el, el + l))

    @staticmethod
    def lcs(S, T):
        m = len(S)
        n = len(T)
        counter = [[0] * (n + 1) for x in range(m + 1)]
        longest = 0
        lcs_set = set()
        for i in range(m):
            for j in range(n):
                if S[i] == T[j]:
                    c = counter[i][j] + 1
                    counter[i + 1][j + 1] = c
                    if c > longest:
                        lcs_set = set()
                        longest = c
                        lcs_set.add(S[i - c + 1:i + 1])
                    elif c == longest:
                        lcs_set.add(S[i - c + 1:i + 1])

        return lcs_set


class Convert:

    @staticmethod
    def from_benchling(data):

        def clean_data(dic):
            for key in dic:
                if isinstance(dic[key], str):
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