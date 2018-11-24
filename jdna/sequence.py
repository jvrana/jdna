"""
Represent linear or circularized nucleotides
"""

import itertools
from collections import defaultdict
from copy import copy, deepcopy
from enum import IntFlag
import textwrap


from jdna.linked_list import Node, DoubleLinkedList, LinkedListMatch
from jdna.utils import random_color
from jdna.alphabet import DNA, UnambiguousDNA, AmbiguousDNA


class SequenceFlags(IntFlag):
    FORWARD = 1
    REVERSE = -1
    TOP = 1
    BOTTOM = -1


class Feature(object):

    def __init__(self, name, type='misc feature', strand=SequenceFlags.FORWARD, color=None):
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


class BindPos(LinkedListMatch):

    def __init__(self, template_bounds, query_bounds, template, query, direction, strand=SequenceFlags.TOP):
        """
        Makes a sequence binding position.

        :param template_bounds_list: list of 2 len tuples containing starts and ends from a template
        :type template_bounds_list: template DoubleLinkedList
        :param query_bounds_list: list of 2 len tuples containing starts and ends from a query
        :type query_bounds_list: query DoubleLinkedList
        :param template: the template
        :type template: DoubleLinkedList
        :param query: the query
        :type query: DoubleLinkedList
        :param direction: If SequenceFlags.FORWARD, the binding position indicates binding forward, to the bottom strand
                            of a dsDNA sequence.
        :type direction: int
        :param strand: If SequenceFlags.BOTTOM, then the query is assumed to be the reverse_complement of the original
                        query
        :type strand: int
        """
        super().__init__(template_bounds, query_bounds, template, query)
        self.direction = direction
        self.strand = strand

        self.anneal = query.copy_slice(*self.query_bounds)
        self.five_prime_overhang = query.copy_slice(None, query_bounds[0].prev())
        self.three_prime_overhang = query.copy_slice(query_bounds[1].next(), None)

        # self.anneal = self.primer[query_span[0]:query_span[1]+1]
        # self.five_prime_overhang = self.primer[:query_span[0]]
        # self.three_prime_overhang = self.primer[query_span[1]+1:]

    # def innitialize(self):
    #     if self.direction == SequenceFlags.REVERSE:
    #         if self.anneal:
    #             self.anneal.reverse_complement()
    #         if self.five_prime_overhang:
    #             self.five_prime_overhang.reverse_complement()
    #         if self.three_prime_overhang:
    #             self.three_prime_overhang.reverse_complement()
    #             length = len(self.three_prime_overhang)
    #         else:
    #             length = 0
    #         self.three_prime_overhang, self.five_prime_overhang = self.five_prime_overhang, self.three_prime_overhang
    #         self.query_span = (self.query_span[0] + length, self.query_span[1] + length)

    @classmethod
    def from_match(cls, linked_list_match, template, query, direction, strand=SequenceFlags.TOP):
        """
        Return a binding pos
        :param linked_list_match: the linked list match
        :type linked_list_match: LinkedListMatch
        :return:
        :rtype:
        """
        return cls(linked_list_match.template_bounds, linked_list_match.query_bounds,
            template, query,
            direction,
            strand=strand
        )

    @property
    def template_anneal(self):
        if self.strand == SequenceFlags.FORWARD:
            return Sequence.new_slice(self.start, self.end)
        else:
            return Sequence.new_slice(self.start, self.end)

    @property
    def query_anneal(self):
        if self.direction == SequenceFlags.FORWARD:
            return Sequence.new_slice(self.query_start, self.query_end)
        else:
            return Sequence.new_slice(self.query_end, self.query_start).reverse_complement()

    def __repr__(self):
        return "<{cls} span={span} direction='{direction}' strand='{strand}' 5'='{five}' anneal='{anneal}' 3'='{three}'>".format(
            cls=self.__class__.__name__,
            span=self.span,
            direction=self.direction,
            strand=self.strand,
            five=self.five_prime_overhang.__repr__(),
            three=self.three_prime_overhang.__repr__(),
            anneal=self.anneal.__repr__()
        )


class Nucleotide(Node):

    __slots__ = ['data', '__next', '__prev', '_features']


    def __init__(self, base):
        super(Nucleotide, self).__init__(base)
        self._features = set()

    @classmethod
    def random(cls):
        """Generate a random sequence"""
        return cls(UnambiguousDNA.random())

    @property
    def base(self):
        return self.data

    def equivalent(self, other):
        return self.base.upper() == other.base.upper()

    def complementary(self, other):
        return self.base.upper() == DNA[other.base].upper()

    def to_complement(self):
        self.data = DNA[self.data]

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
    #         n.ffeatures[feature] += delta_i

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
        copied._features = deepcopy(self.features)
        return copied


class Sequence(DoubleLinkedList):

    NODE_CLASS = Nucleotide
    counter = itertools.count()

    def __init__(self, sequence=None, first=None, name=None, description=''):
        super(Sequence, self).__init__(data=sequence, first=first)
        self.name = name
        self.description = description
        self._global_id = next(Sequence.counter)

    @property
    def global_id(self):
        return self._global_id

    @classmethod
    def random(cls, length):
        """Generate a random sequence"""
        seq = ""
        for i in range(length):
            seq += UnambiguousDNA.random().upper()
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
                    if positions[0][0] == 0 and positions[-1][-1] == l - 1:
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

    def create_feature(self, start, end, name, feature_type):
        f = Feature(name, feature_type)
        self.add_feature(start, end, f)
        return f

    def annotate(self, start, end, name, feature_type=None):
        if feature_type is None:
            feature_type = 'misc'
        return self.create_feature(start, end, name, feature_type)

    def complement(self):
        curr = self.head
        visited = set()
        while curr and curr not in visited:
            visited.add(curr)
            curr.to_complement()
            curr = next(curr)
        return self

    def c(self):
        return self.complement()

    def reverse_complement(self):
        self.reverse()
        self.complement()
        return self

    def rc(self):
        return self.reverse_complement()

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

    # def anneal_to_bottom_strand(self, other, min_bases=10):
    #     for match in self.find_iter(other,
    #                                 min_query_length=min_bases,
    #                                 direction=self.Direction.REVERSE, ):
    #         yield match
    #
    # def anneal_to_top_strand(self, other, min_bases=10):
    #     for match in self.find_iter(other,
    #                                 min_query_length=min_bases,
    #                                 protocol=lambda x, y: x.complementary(y)):
    #         yield match

    def anneal_forward(self, other, min_bases=10):
        for match in self.find_iter(other, min_query_length=min_bases,
                                    direction=self.Direction.REVERSE):
            yield BindPos.from_match(match, self, other, direction=self.Direction.FORWARD)

    def anneal_reverse(self, other, min_bases=10):
        for match in self.find_iter(other,
                                    min_query_length=min_bases,
                                    direction=(1, -1),
                                    protocol=lambda x, y: x.complementary(y)
                                    ):
            yield BindPos.from_match(match, self, other, direction=self.Direction.REVERSE)

    def anneal(self, ssDNA, min_bases=10):
        """Simulate annealing a single stranded piece of DNA to a double_stranded template"""
        for match in self.anneal_forward(ssDNA, min_bases=min_bases):
            yield match
        for match in self.anneal_reverse(ssDNA, min_bases=min_bases):
            yield match

    def dsanneal(self, dsDNA, min_bases=10):
        """Simulate annealing a double stranded piece of DNA to a double_stranded template"""
        for binding in self.anneal(dsDNA, min_bases=min_bases):
            yield binding
        for binding in self.anneal(dsDNA.copy().reverse_complement(), min_bases=min_bases):
            binding.strand = SequenceFlags.BOTTOM
            yield binding

    def __repr__(self):
        max_width = 30
        replace = '...'
        display = int((max_width - len(replace))/2.0)
        s = str(self)
        if len(s) > display*2:
            # diff = display*2 - len(s)
            s = s[:display] + '...' + s[-display:]
        return "Sequence('{}')".format(s)