"""
Represent linear or circularized nucleotides
"""

from copy import copy, deepcopy
import itertools
from collections import defaultdict
from jdna.linked_list import Node, DoubleLinkedList
from jdna.utils import random_color
from enum import Enum
import random

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

    def __copy__(self):
        return self.__class__(self.name, self.type, self.strand, self.color)

    @property
    def nodes(self):
        return self._nodes

    def segments(self):
        return Sequence.segments(self.nodes)
        # visited = set()
        # pairs = set()
        # stop = lambda x: x not in self._nodes
        # for n in self._nodes:
        #     if n not in visited:
        #         tail = n
        #         for tail in n.fwd(stop_criteria=stop):
        #             visited.add(tail)
        #         head = n
        #         for head in n.rev(stop_criteria=stop):
        #             visited.add(head)
        #         pairs.add((head, tail))
        # return pairs

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

    BASES = dict(list(zip(
        ['a', 't', 'c', 'g', 'A', 'T', 'C', 'G', 'n', 'N'],
        ['t', 'a', 'g', 'c', 'T', 'A', 'G', 'C', 'n', 'N']
    )))

    def __init__(self, base):
        super(Nucleotide, self).__init__(base)
        self._features = set()

    @classmethod
    def random(cls):
        """Generate a random sequence"""
        return cls(random.choice(cls.BASES))

    def base(self):
        return self.data

    def equivalent(self, other):
        return str(self.base()).upper() == str(other.base()).upper()

    def to_complement(self):
        self.data = self.BASES[self.data]

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
        if cut_prev:
            nxt = super().cut_prev()
        else:
            nxt = super().cut_next()
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
        if n1 is None or n2 is None:
            return
        if fuse_condition is None:
            fuse_condition = cls._default_fuse_condition

        if not (n1.next() is n2 and n2.prev() is n1):
            raise Exception("Cannot fuse non-consecutive features")

        for f1 in set(n1.features):
            for f2 in set(n2.features):
                if f1 is not f2 and fuse_condition(f1, f2):
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
    counter = itertools.count()

    def __init__(self, sequence=None, first=None, name=None, description=''):
        super(Sequence, self).__init__(data=sequence, first=first)
        self.name = name
        self.description = description
        self._global_id = next(Sequence.counter)

    @classmethod
    def random(cls, length):
        """Generate a random sequence"""
        seq = ""
        for i in range(length):
            seq += random.choice(list(cls.NODE_CLASS.BASES.keys()))
        if seq == '':
            return None
        return cls(sequence=seq)

    @property
    def features(self):
        features_set = set()
        for i, n in enumerate(self):
            features_set.update(n.features)
        return features_set

    def feature_positions(self, with_nodes=False):
        index = 0
        feature_pos = defaultdict(list)
        feature_nodes = defaultdict(list)
        l = len(self)
        for n in self:
            for f in n.features:
                if feature_pos[f] and feature_pos[f][-1][-1] + 1 == index:
                    feature_pos[f][-1][-1] = index
                    feature_nodes[f][-1][-1] = n
                else:
                    feature_pos[f].append([index, index])
                    feature_nodes[f].append([n, n])
            index += 1

        # capture features that span the origin
        if self.cyclic:
            for k in feature_pos:
                positions = feature_pos[k]
                nodes = feature_nodes[k]
                if len(nodes) > 1:
                    if positions[0][0] == 0 and positions[-1][-1] == l-1:
                        nodes[0][0] = nodes[-1][0]
                        positions[0][0] = positions[-1][0]
                        nodes.pop()
                        positions.pop()

        if with_nodes:
            return feature_pos, feature_nodes
        return feature_pos

    def feature_nodes(self):
        return self.feature_positions(with_nodes=True)[-1]

    def add_feature(self, i, j, feature):
        feature_nts = list(self.inclusive_range(i, j))
        if feature_nts[-1] is not self[j]:
            if not self.cyclic:
                raise IndexError("Cannot add feature to {} to linear dna with bounds {}".format(
                    (i, j),
                    (0, len(self))
                ))
            else:
                raise IndexError("Cannot add feature to {}".format(
                    (i, j)
                ))
        for n in self.inclusive_range(i, j):
            n.add_feature(feature)
        return feature

    def add_multipart_feature(self, positions, feature):
        for i, j in positions:
            self.add_feature(i, j, feature)
        return feature

    def print_features(self):
        raise NotImplementedError()

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

    def __copy__(self):
        feature_positions = self.feature_positions()
        copied = super(Sequence, self).__copy__()

        for feature, positions in feature_positions.items():
            copied.add_multipart_feature(positions, copy(feature))
        return copied

    # def __repr__(self):
    #     return str(self)

