from copy import copy, deepcopy
import itertools
from collections import defaultdict
from jdna.linked_list import Node, DoubleLinkedList
from jdna.utils import random_color
from enum import Enum

class STRAND(Enum):

    FORWARD = 1
    REVERSE = -1

class Feature(object):

    def __init__(self, name, type='misc feature', strand=STRAND.FORWARD, color=None):
        self.name = name
        self.type = type
        self.strand = strand
        if color is None:
            color = random_color()
        self.color = color
        self._nodes = set()

    def __str__(self):
        return '{} {}'.format(self.name, self.type)

    def __repr__(self):
        return str(self)

    @property
    def nodes(self):
        return self._nodes

    @property
    def segments(self):
        visited = set()
        pairs = set()
        stop = lambda x: x not in self._nodes
        for n in self._nodes:
            if n not in visited:
                tail = n
                for tail in n.fwd(stop_criteria=stop):
                    visited.add(tail)
                head = n
                for head in n.rev(stop_criteria=stop):
                    visited.add(head)
                pairs.add((head, tail))
        return pairs

    def is_multipart(self):
        if len(self.segments) > 1:
            return True
        return False

    def _bind(self, nodes):
        for n in nodes:
            self._nodes.add(n)

    def _unbind(self, nodes):
        for n in nodes:
            if n in self._nodes:
                self._nodes.remove(n)


class Nucleotide(Node):

    base_pairing = dict(list(zip(
        ['a', 't', 'c', 'g', 'A', 'T', 'C', 'G', 'n', 'N'],
        ['t', 'a', 'g', 'c', 'T', 'A', 'G', 'C', 'n', 'N']
    )))

    def __init__(self, base):
        super(Nucleotide, self).__init__(base)
        self._features = set()

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

    @property
    def features(self):
        return self._features

    def add_feature(self, feature):
        self.features.add(feature)
        feature._bind([self])

    def remove_feature(self, feature):
        self.features.remove(feature)
        feature._unbind([self])

    def feature_fwd(self, feature):
        stop = lambda x: feature not in x.features
        return self._propogate(lambda x: x.next(), stop_criteria=stop)

    def feature_rev(self, feature):
        stop = lambda x: feature not in x.features
        return self._propogate(lambda x: x.prev(), stop_criteria=stop)

    def replace_feature(self, old_feature, new_feature):
        self.features[new_feature] = self.features[old_feature]
        self.remove_feature(old_feature)

    def copy_features_from(self, other):
        for f in other.features:
            i = other.features[f]
            if f not in self.features:
                self.add_feature(f, i)
        self._remove_overlapping_features()

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
    def _default_fuse_condition(f1, f2):
        return f1.name == f2.name

    @classmethod
    def fuse_features(cls, n1, n2, fuse_condition=None):

        if fuse_condition is None:
            fuse_condition = cls._default_fuse_condition

        if not (n1.next is n2 and n2.prev is n1):
            n1_next = id(n1.next)
            n2_prev = id(n2.prev)
            n1_id = id(n1)
            n2_id = id(n2)
            n = n1.next
            _n = n2.prev
            raise Exception("Cannot fuse non-consecutive features")

        for f1 in n1.features:
            for f2 in n2.features:
                if fuse_condition(f1, f2):
                    for n in n2.feature_fwd(f2):
                        n.add_feature(f1)
                        n.remove_feature(f2)


    #
    # @staticmethod
    # def fuse_features(n1, n2):
    #     if n1 is None:
    #         return
    #     if n2 is None:
    #         return
    #
    #     delset = set()
    #
    #     for f1 in n1.features:
    #         for f2 in n2.features:
    #             f1_pos = n1.features[f1]
    #             f2_pos = n2.features[f2]
    #             f1_copy = copy(f1)
    #             # same name & consecutive position
    #             if f1 is f2:
    #                 continue
    #             if f1.name == f2.name and f1_pos + 1 == f2_pos:
    #                 delset.add((f1, f2, f1_copy))
    #     for f1, f2, f1_copy in delset:
    #         for n in n1.feature_rev(f1):
    #             try:
    #                 n.replace_feature(f1, f1_copy)
    #             except KeyError:
    #                 pass
    #         for n in n2.feature_fwd(f2):
    #             try:
    #                 n.replace_feature(f2, f1_copy)
    #             except KeyError:
    #                 pass

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

    NODE_CLASS = Nucleotide

    def __init__(self, first=None, sequence=None, name='unknown'):
        super(Sequence, self).__init__(first=first, sequence=sequence)
        self.name = name
        self.description = ''

    @property
    def features(self):
        features_set = set()
        for n in self:
            features_set = features_set.union(n.features)
        return features_set

    def feature_positions(self, with_nodes=False):
        index = 0
        feature_dict = defaultdict(list)
        feature_nodes = defaultdict(list)
        for n in self:
            for f in n.features:
                if feature_dict[f] and feature_dict[f][-1][-1] + 1 == index:
                    feature_dict[f][-1][-1] = index
                    feature_nodes[f][-1][-1] = n
                else:
                    feature_dict[f].append([index, index])
                    feature_nodes[f].append([n, n])
            index += 1
        if with_nodes:
            return feature_dict, feature_nodes
        return feature_dict

    def feature_nodes(self):
        return self.feature_positions(with_nodes=True)[-1]

    def add_feature(self, i, j, feature):
        s = self.nodes
        nodes = self.nodes
        l = len(nodes)
        start = s[i]
        end = s[j]
        feature_nts = list(start.fwd(stop_node=end))
        if end != feature_nts[-1]:
            raise IndexError("")
        for i, n in enumerate(feature_nts):
            n.add_feature(feature)
        feature.length = len(feature_nts)
        return feature

    def add_multipart_feature(self, positions, feature):
        for i, j in positions:
            self.add_feature(i, j, feature)
        return feature

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

    def find_feature_by_name(self, name):
        found = []
        for feature in self.feature_positions():
            if feature.name == name:
                found.append(feature)
        return found

    def create_feature(self, name, feature_type, start, end):
        f = Feature(name, feature_type)
        self.add_feature(start, end, f)
        return f

    def complement(self):
        curr = self.head
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
        fragments = [Sequence(first=f.head) for f in fragments]
        return fragments

    def chop_off_fiveprime(self, i):
        if self.cyclic:
            raise IndexError('Cannot chop a cyclic sequence.')
        return self.cut(i)[-1]

    def chop_off_threeprime(self, i):
        if self.cyclic:
            raise IndexError('Cannot chop a cyclic sequence.')
        return self.cut(i, cut_prev=False)[0]

    def fuse(self, seq):
        f = self.head.find_last()
        l = seq.head
        f.set_next(l)
        return self

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        copied = super(Sequence, self).__copy__()
        for feature, positions in self.feature_positions().items():
            self.add_multipart_feature(positions, feature)
        return copied

    def __add__(self, other):
        return copy(self).fuse(copy(other))

    # def __repr__(self):
    #     return str(self)

