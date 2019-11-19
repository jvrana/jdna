"""Represent linear or circularized nucleotides."""
import itertools
import re
from collections import defaultdict
from copy import copy
from typing import Any
from typing import Callable
from typing import Dict
from typing import List
from typing import Sequence as TypingSequence
from typing import Tuple
from typing import Union

import primer3
from Bio import Restriction

from jdna.align import AlignInterface
from jdna.alphabet import AmbiguousDNA
from jdna.alphabet import UnambiguousDNA
from jdna.format import format_sequence
from jdna.io import IOInterface
from jdna.linked_list import DoubleLinkedList
from jdna.linked_list import LinkedListMatch
from jdna.linked_list import Node
from jdna.utils import random_color
from jdna.viewer import SequenceViewer
from jdna.viewer import StringColumn
from jdna.viewer import ViewerAnnotationFlag


class SequenceFlags:
    """Constants/Flags for sequences."""

    FORWARD = 1
    REVERSE = -1
    TOP = 1
    BOTTOM = -1


class Feature:
    """An annotation for a sequence."""

    def __init__(self, name, type=None, strand=None, color=None):
        self.name = name
        if type is None:
            type = "misc"
        self.type = type
        if strand is None:
            strand = SequenceFlags.FORWARD
        self.strand = strand
        if color is None:
            color = random_color()
        self.color = color
        # self._nodes = set()

    def reverse(self) -> "Sequence":
        self.strand = -1 * self.strand
        return self

    def __str__(self) -> str:
        return "<Feature name='{name}' type='{tp}' color='{color}'".format(
            name=self.name, tp=self.type, color=self.color
        )

    def __repr__(self) -> str:
        return str(self)

    def __copy__(self) -> "Sequence":
        return self.__class__(self.name, self.type, self.strand, self.color)

    # @property
    # def nodes(self):
    #     return self._nodes
    #
    # def segments(self):
    #     return Sequence.segments(self.nodes)
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

    def is_multipart(self) -> bool:
        if len(self.segments) > 1:
            return True
        return False

    def _bind(self, nodes):
        pass
        # for n in nodes:
        #     self._nodes.add(n)

    # def _unbind(self, nodes):
    #     for n in nodes:
    #         if n in self._nodes:
    #             self._nodes.remove(n)


class BindPos(LinkedListMatch):
    def __init__(
        self,
        template_bounds: Tuple["Nucleotide", "Nucleotide"],
        query_bounds: Tuple["Nucleotide", "Nucleotide"],
        template: "Sequence",
        query: "Sequence",
        direction: int,
        strand=SequenceFlags.TOP,
    ):
        """Makes a sequence binding position.

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

        if self.direction == SequenceFlags.REVERSE:
            self.anneal = query.copy_slice(*self.query_bounds[::-1])
            self.five_prime_overhang = query.new_slice(None, self.query_end.prev())
            self.three_prime_overhang = query.new_slice(self.query_start.next(), None)
        else:
            self.anneal = query.copy_slice(*self.query_bounds)
            self.five_prime_overhang = query.new_slice(None, self.query_start.prev())
            self.three_prime_overhang = query.new_slice(self.query_end.next(), None)
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
    def from_match(
        cls, linked_list_match, template, query, direction, strand=SequenceFlags.TOP
    ):
        """Return a binding pos.

        :param linked_list_match: the linked list match
        :type linked_list_match: LinkedListMatch
        :return:
        :rtype:
        """
        return cls(
            linked_list_match.template_bounds,
            linked_list_match.query_bounds,
            template,
            query,
            direction,
            strand=strand,
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
            return Sequence.new_slice(
                self.query_end, self.query_start
            ).reverse_complement()

    def __repr__(self):
        return "<{cls} span={span} direction='{direction}' strand='{strand}' 5'='{five}' anneal='{anneal}' 3'='{three}'>".format(
            cls=self.__class__.__name__,
            span=self.span,
            direction=self.direction,
            strand=self.strand,
            five=self.five_prime_overhang.__repr__(),
            three=self.three_prime_overhang.__repr__(),
            anneal=self.anneal.__repr__(),
        )


class Nucleotide(Node):
    """Represents a biological nucleotide.

    Serves a :class:`Node` in teh :class:`Sequence` object.
    """

    __slots__ = ["data", "__next", "__prev", "_features", "alphabet"]

    def __init__(self, base, alphabet=AmbiguousDNA):
        """Nucleotide constructor.

        :param base: base as a single character string
        :type base: basestring
        """
        super().__init__(base)
        self._features = set()
        self.alphabet = alphabet

    @classmethod
    def random(cls):
        """Generate a random sequence."""
        return cls(UnambiguousDNA.random())

    @property
    def base(self):
        return self.data

    def equivalent(self, other) -> bool:
        return self.alphabet.compare(self.base, other.base)

    def complementary(self, other) -> bool:
        return self.base.upper() == AmbiguousDNA[other.base].upper()

    def to_complement(self):
        self.data = AmbiguousDNA[self.data]

    def set_next(self, nucleotide):
        self.cut_next()
        super().set_next(nucleotide)
        Nucleotide.fuse_features(self, nucleotide)

    def set_prev(self, nucleotide):
        self.cut_prev()
        super().set_prev(nucleotide)
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
        return feature

    def remove_feature(self, feature):
        self.features.remove(feature)

    def feature_fwd(self, feature):
        def stop(x):
            return feature not in x.features

        return self._propogate(lambda x: x.next(), stop_criteria=stop)

    def feature_rev(self, feature):
        def stop(x):
            return feature not in x.features

        return self._propogate(lambda x: x.prev(), stop_criteria=stop)

    def replace_feature(self, old_feature, new_feature):
        self.features[new_feature] = self.features[old_feature]
        self.remove_feature(old_feature)

    def copy_features_from(self, other):
        for f in other.features:
            if f not in self.features:
                self.add_feature(f)
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

    def _clear_features(self):
        self._features = set()

    def copy(self):
        copied = super().copy()
        copied._features = set()
        for f in self.features:
            copied.add_feature(f)
        return copied


class Sequence(DoubleLinkedList):
    """Represents a biological sequence as a double linked list.

    Can be annotated with features.
    """

    class DEFAULTS:
        """Sequence defaults."""

        MIN_ANNEAL_BASES = 13
        FOREGROUND_COLORS = ["blue", "red"]
        BACKGROUND_COLORS = None
        ALPHABET = AmbiguousDNA

    FORWARD = SequenceFlags.FORWARD
    REVERSE = SequenceFlags.REVERSE
    TOP = SequenceFlags.TOP
    BOTTOM = SequenceFlags.BOTTOM
    NODE_CLASS = Nucleotide
    counter = itertools.count()

    def __init__(
        self,
        sequence: TypingSequence[Any] = None,
        first: Nucleotide = None,
        name: str = None,
        description: str = "",
        metadata: dict = None,
        cyclic: bool = False,
        alphabet=DEFAULTS.ALPHABET,
    ):
        """

        :param sequence: sequence string
        :type sequence: basestring
        :param first: optional first Nucleotide to use as the 'head' to this Sequence
        :type first: Nucleotide
        :param name: optional name of the sequence
        :type name: basestring
        :param description: optional description of the sequence
        :type description: basestring
        :param metadata: additional sequence metadata
        :type metadata: dict
        :param cyclic: whether to make the sequence circular
        :type cyclic: bool
        :param alphabet: the base pair alphabet of this sequence which used for complementary and comparisons
                        (default: AmbiguousDNA)
        :type alphabet: jdna.alphabet.Alphabet
        """

        self.alphabet = alphabet
        super().__init__(data=sequence, first=first, cyclic=cyclic)
        if name is None:
            name = ""
        self.name = name
        self.description = description
        if metadata is None:
            metadata = dict()
        self.metadata = metadata
        self._global_id = next(Sequence.counter)
        self._io = self.IO.instance(self)
        self._align = self.Align.instance(self)
        if cyclic:
            self.cyclic = cyclic

    def new_node(self, data):
        return self.NODE_CLASS(data, alphabet=self.alphabet)

    @property
    def io(self):
        return self._io

    @property
    def align(self):
        return self._align

    @property
    def global_id(self):
        return self._global_id

    @classmethod
    def random(cls, length):
        """Generate a random sequence."""
        seq = ""
        for i in range(length):
            seq += UnambiguousDNA.random().upper()
        if seq == "":
            return cls.empty()
        return cls(sequence=seq)

    @property
    def features_list(self):
        """Returns set of features contained in sequence.

        :return: set of features in this sequence
        :rtype: set
        """
        features_set = set()
        for i, n in enumerate(self):
            features_set.update(n.features)
        return tuple(features_set)

    @property
    def features(self, with_nodes=False):
        """Return a list of feature positions.

        :param with_nodes: if True, will return a tuple composed of a feature to
                            position dictionary and a feature to
                            start and end node. If False, will just return a feature to
                            position dictionary
        :type with_nodes: bool
        :return: feature positions dictionary OR tuple of feature positions dictionary
                            and feature node dictionary
        :rtype: tuple
        """
        index = 0
        feature_pos = defaultdict(list)
        feature_nodes = defaultdict(list)
        length = len(self)
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
                    if positions[0][0] == 0 and positions[-1][-1] == length - 1:
                        nodes[0][0] = nodes[-1][0]
                        positions[0][0] = positions[-1][0]
                        nodes.pop()
                        positions.pop()

        if with_nodes:
            return feature_pos, feature_nodes
        return feature_pos

    def feature_nodes(self):
        return self.features(with_nodes=True)[-1]

    def add_feature(self, start, end, feature):
        """Add a feature to the start and end positions (inclusive)

        :param start: start
        :type start: int
        :param end: end (inclusive)
        :type end: int
        :param feature: the feature to add
        :type feature: Feature
        :return: the added feature
        :rtype: Feature
        """
        feature_nts = list(self.inclusive_range(start, end))
        if end and feature_nts[-1] is not self[end]:
            if not self.cyclic:
                raise IndexError(
                    "Cannot add feature to {} to linear dna with bounds {}".format(
                        (start, end), (0, len(self))
                    )
                )
            else:
                raise IndexError("Cannot add feature to {}".format((start, end)))
        for n in self.inclusive_range(start, end):
            n.add_feature(feature)
        return feature

    def add_multipart_feature(self, positions, feature):
        """Add a multi-part feature (i.e. a disjointed feature)

        :param positions: list of start and ends as tuples ([(1,100), (110,200)]
        :type positions: list
        :param feature: the feature to add
        :type feature: Feature
        :return: the added feature
        :rtype: Feature
        """
        for i, j in positions:
            self.add_feature(i, j, feature)
        return feature

    # def print_features(self):
    #     raise NotImplementedError()

    def find_feature_by_name(self, name):
        """Find features by name.

        :param name: feature name
        :type name: basestring
        :return: list of features
        :rtype: list
        """
        found = []
        for feature in self.features:
            if feature.name == name:
                found.append(feature)
        return found

    def annotate(self, start, end, name, feature_type=None, color=None, strand=None):
        """Annotate a regions.

        :param start: start
        :type start: int
        :param end: end (inclusive)
        :type end: end
        :param name: feature name
        :type name: basestring
        :param feature_type: feature type (default=misc)
        :type feature_type: basestring
        :param color: optional feature color
        :type color: basestring
        :return: new feature
        :rtype: Feature
        """
        return self.add_feature(
            start, end, Feature(name, feature_type, strand=strand, color=color)
        )

    def reverse(self):
        features_set = set()
        if self.is_empty():
            return self
        nodes = self.nodes
        for s in nodes:
            s.swap()
            features_set.update(s.features)
        for f in features_set:
            f.reverse()
        self.head = nodes[-1]
        return self

    def complement(self):
        """Complement the sequence in place."""
        if self.is_empty():
            return self
        curr = self.head
        visited = set()
        while curr and curr not in visited:
            visited.add(curr)
            curr.to_complement()
            curr = next(curr)
        return self

    def c(self):
        """Complement the sequence in place."""
        return self.complement()

    def reverse_complement(self):
        """Reverse complement the sequence in place."""
        self.reverse()
        self.complement()
        return self

    def rc(self):
        """Reverse complement the sequence in place."""
        return self.reverse_complement()

    def cut(self, i, cut_prev=True):
        fragments = super().cut(i, cut_prev)
        fragments = [Sequence(first=f.head) for f in fragments]
        return fragments

    def clear_features(self):
        for n in self:
            n._clear_features()

    def __copy__(self):
        copied = super().__copy__()
        copied.name = self.name
        copied._global_id = next(self.counter)
        copied.clear_features()
        feature_positions = self.features
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

    def anneal_forward(self, other, min_bases=DEFAULTS.MIN_ANNEAL_BASES, depth=None):
        """Anneal a sequence in the forward direction."""
        for match in self.find_iter(
            other,
            min_query_length=min_bases,
            direction=self.Direction.REVERSE,
            depth=depth,
        ):
            yield BindPos.from_match(
                match, self, other, direction=self.Direction.FORWARD
            )

    def anneal_reverse(self, other, min_bases=DEFAULTS.MIN_ANNEAL_BASES, depth=None):
        """Anneal a sequence in the reverse direction."""
        for match in self.find_iter(
            other,
            min_query_length=min_bases,
            direction=(1, -1),
            protocol=lambda x, y: x.complementary(y),
            depth=depth,
        ):
            yield BindPos.from_match(
                match, self, other, direction=self.Direction.REVERSE
            )

    def anneal(self, ssDNA, min_bases=DEFAULTS.MIN_ANNEAL_BASES, depth=None):
        """Simulate annealing a single stranded piece of DNA to a
        double_stranded template."""
        for match in self.anneal_forward(ssDNA, min_bases=min_bases, depth=depth):
            yield match
        for match in self.anneal_reverse(ssDNA, min_bases=min_bases, depth=depth):
            yield match

    def dsanneal(self, dsDNA, min_bases=DEFAULTS.MIN_ANNEAL_BASES, depth=None):
        """Simulate annealing a double stranded piece of DNA to a
        double_stranded template."""
        for binding in self.anneal(dsDNA, min_bases=min_bases, depth=depth):
            yield binding
        for binding in self.anneal(
            dsDNA.copy().reverse_complement(), min_bases=min_bases, depth=depth
        ):
            binding.strand = SequenceFlags.BOTTOM
            yield binding

    def format(self, width=75, spacer=""):
        return format_sequence(str(self), width=width, spacer=spacer)

    @classmethod
    def _apply_features_to_view(cls, sequence, view):
        for feature, positions in sequence.features.items():
            for pos in positions:
                direction = None
                if feature.strand == SequenceFlags.FORWARD:
                    direction = ViewerAnnotationFlag.FORWARD
                elif feature.strand == SequenceFlags.REVERSE:
                    direction = ViewerAnnotationFlag.REVERSE
                view.annotate(
                    pos[0],
                    pos[1],
                    label=feature.name,
                    fill=direction,
                    background=feature.color,
                )

    def view_bindings(self, bindings, view=None):
        if view is None:
            view = self.view(complement=True)
        for b in bindings:
            anneal = b.anneal
            primer_sequence = b.five_prime_overhang + anneal + b.three_prime_overhang
            annotation = StringColumn(
                [
                    str(primer_sequence),
                    " " * len(b.five_prime_overhang)
                    + "|" * len(anneal)
                    + " " * len(b.three_prime_overhang),
                ]
            )
            if b.direction == Sequence.FORWARD:
                view.annotate(b.span[0], b.span[1], annotation)
            if b.direction == Sequence.REVERSE:
                view.annotate(b.span[0], b.span[1], annotation.flip()[::-1], top=False)
        return view

    def view(
        self,
        indent=10,
        width=85,
        spacer=None,
        complement=False,
        features=True,
        **kwargs
    ):
        """Create a :class:`SequenceViewer` instance from this sequence.
        Printing the view object with annotations and complement will produce
        an output similar to the following:

        .. code::


            > "Unnamed" (550bp)


                                                                        ----------------GFP----------------
                                                                        |<START
                                                                        ----      -----------RFP-----------
            0         CCCAGGACTAGCGACTTTCCGTAACGCGACCTAACACCGGCCGTTCCTTCGAGCCAGGCAAATGTTACGTCACTTCCTTAGATTT
                      GGGTCCTGATCGCTGAAAGGCATTGCGCTGGATTGTGGCCGGCAAGGAAGCTCGGTCCGTTTACAATGCAGTGAAGGAATCTAAA

                      ------GFP------
                      -----------------------------------------RFP-----------------------------------------
            85        TGAACAGCGCCGTACCCCGATATGATATTTAGATATATAGCAGTTACACTTGGGGTTGCTATGGACTTAGATCTGCTGTATGTTT
                      ACTTGTCGCGGCATGGGGCTATACTATAAATCTATATATCGTCAATGTGAACCCCAACGATACCTGAATCTAGACGACATACAAA

                      -----------------------------------------RFP-----------------------------------------
            170       TCTTACCTTCCGCATCAGGGGACAATTCGCCAGTAGAATTCAGTTTGTGCGTGAGAACATAAGATTGAATCCCACGCAGGCACAA
                      AGAATGGAAGGCGTAGTCCCCTGTTAAGCGGTCATCTTAAGTCAAACACGCACTCTTGTATTCTAACTTAGGGTGCGTCCGTGTT

                      ---------------------RFP----------------------
            255       GCAGGGCGGGCAGACTCTATAGGTCCTAAGACCCTGAGACTGCGTCCTCAAGATACAGGTTAACAATCCCCGTATGGAGCCGTTC
                      CGTCCCGCCCGTCTGAGATATCCAGGATTCTGGGACTCTGACGCAGGAGTTCTATGTCCAATTGTTAGGGGCATACCTCGGCAAG

            340       TTAGCATGACCCGACAGGTGGGCTTGGCTCGCGTAAGTTGAGTGTTGCAGATACCTGCTGCTGCGCGGTCTAGGGGGAATCGCCG
                      AATCGTACTGGGCTGTCCACCCGAACCGAGCGCATTCAACTCACAACGTCTATGGACGACGACGCGCCAGATCCCCCTTAGCGGC

            425       ATTTTGACGTAGGATCGGTAATGGGCAGTAAACCCGCAACTATTTTCAGCACCAGATGCAAGTTTCCCTAGAAAGCGTCATGGTT
                      TAAAACTGCATCCTAGCCATTACCCGTCATTTGGGCGTTGATAAAAGTCGTGGTCTACGTTCAAAGGGATCTTTCGCAGTACCAA

            510       TGCAATCTCCTTAGGTCACAGCAAACATAGCAGCCCCTGT
                      ACGTTAGAGGAATCCAGTGTCGTTTGTATCGTCGGGGACA

        :param indent: indent between left column and base pairs view windo
        :type indent: int
        :param width: width of the view window
        :type width: int
        :param spacer: string to intersperse between sequence rows (default is newline)
        :type spacer: basestring
        :param complement: whether to include the complementary strand in the view
        :type complement: bool
        :param features: whether to include annotations/features in the view instance
        :type features: bool
        :return: the viewer object
        :rtype: SequenceViewer
        """
        if indent is None:
            indent = 10

        if width is None:
            width = 85

        seqs = [self]
        colors = self.DEFAULTS.FOREGROUND_COLORS[0]
        if complement:
            seqs.append(self.copy().complement())
            colors = self.DEFAULTS.FOREGROUND_COLORS
        if spacer is None:
            if complement:
                spacer = "\n"
            else:
                spacer = ""
        viewer = SequenceViewer(
            seqs,
            name=self.name,
            description=self.description,
            indent=indent,
            width=width,
            spacer=spacer,
            foreground_colors=colors,
            **kwargs
        )
        viewer.metadata.update(self.metadata)
        if features:
            self._apply_features_to_view(self, viewer)
        return viewer

    def upper(self):
        copied = self
        for n in copied:
            n.data = n.data.upper()
        return copied

    def lower(self):
        copied = self
        for n in copied:
            n.data = n.data.lower()
        return copied

    def print(
        self,
        indent=None,
        width=None,
        spacer=None,
        complement=False,
        features=True,
        **kwargs
    ):
        """Create and print a :class:`SequenceViewer` instance from this
        sequence. Printing the view object with annotations and complement will
        produce an output similar to the following:

        .. code::


            > "Unnamed" (550bp)


                                                                        ----------------GFP----------------
                                                                        |<START
                                                                        ----      -----------RFP-----------
            0         CCCAGGACTAGCGACTTTCCGTAACGCGACCTAACACCGGCCGTTCCTTCGAGCCAGGCAAATGTTACGTCACTTCCTTAGATTT
                      GGGTCCTGATCGCTGAAAGGCATTGCGCTGGATTGTGGCCGGCAAGGAAGCTCGGTCCGTTTACAATGCAGTGAAGGAATCTAAA

                      ------GFP------
                      -----------------------------------------RFP-----------------------------------------
            85        TGAACAGCGCCGTACCCCGATATGATATTTAGATATATAGCAGTTACACTTGGGGTTGCTATGGACTTAGATCTGCTGTATGTTT
                      ACTTGTCGCGGCATGGGGCTATACTATAAATCTATATATCGTCAATGTGAACCCCAACGATACCTGAATCTAGACGACATACAAA

                      -----------------------------------------RFP-----------------------------------------
            170       TCTTACCTTCCGCATCAGGGGACAATTCGCCAGTAGAATTCAGTTTGTGCGTGAGAACATAAGATTGAATCCCACGCAGGCACAA
                      AGAATGGAAGGCGTAGTCCCCTGTTAAGCGGTCATCTTAAGTCAAACACGCACTCTTGTATTCTAACTTAGGGTGCGTCCGTGTT

                      ---------------------RFP----------------------
            255       GCAGGGCGGGCAGACTCTATAGGTCCTAAGACCCTGAGACTGCGTCCTCAAGATACAGGTTAACAATCCCCGTATGGAGCCGTTC
                      CGTCCCGCCCGTCTGAGATATCCAGGATTCTGGGACTCTGACGCAGGAGTTCTATGTCCAATTGTTAGGGGCATACCTCGGCAAG

            340       TTAGCATGACCCGACAGGTGGGCTTGGCTCGCGTAAGTTGAGTGTTGCAGATACCTGCTGCTGCGCGGTCTAGGGGGAATCGCCG
                      AATCGTACTGGGCTGTCCACCCGAACCGAGCGCATTCAACTCACAACGTCTATGGACGACGACGCGCCAGATCCCCCTTAGCGGC

            425       ATTTTGACGTAGGATCGGTAATGGGCAGTAAACCCGCAACTATTTTCAGCACCAGATGCAAGTTTCCCTAGAAAGCGTCATGGTT
                      TAAAACTGCATCCTAGCCATTACCCGTCATTTGGGCGTTGATAAAAGTCGTGGTCTACGTTCAAAGGGATCTTTCGCAGTACCAA

            510       TGCAATCTCCTTAGGTCACAGCAAACATAGCAGCCCCTGT
                      ACGTTAGAGGAATCCAGTGTCGTTTGTATCGTCGGGGACA

        :param indent: indent between left column and base pairs view windo
        :type indent: int
        :param width: width of the view window
        :type width: int
        :param spacer: string to intersperse between sequence rows (default is newline)
        :type spacer: basestring
        :param complement: whether to include the complementary strand in the view
        :type complement: bool
        :param include_annotations: whether to include annotations/features in the view instance
        :type include_annotations: bool
        :return: the viewer object
        :rtype: SequenceViewer
        """
        self.view(
            indent=indent,
            width=width,
            spacer=spacer,
            complement=complement,
            features=features,
            **kwargs
        ).print()

    def tm(self):
        """Calculate the Tm of this sequence using primer3 defaults.

        :return: the tm of the sequence
        :rtype: float
        """
        return primer3.calcTm(str(self).upper())

    def json(self):
        """Print sequence to a json dictionary."""
        annotations = []
        for feature, positions in self.features.items():
            for start, end in positions:
                annotations.append(
                    {
                        "start": start,
                        "end": end + 1,
                        "name": feature.name,
                        "color": feature.color,
                        "type": feature.type,
                        "strand": feature.strand,
                    }
                )

        return {
            "name": self.name,
            "isCircular": self.cyclic,
            "length": len(self),
            "bases": str(self),
            "annotations": annotations,
        }

    @classmethod
    def load(cls, data):
        """Load a sequence from a json formatted dictionary."""
        sequence = cls(data["bases"], name=data["name"])
        sequence.cyclic = data["isCircular"]
        sequence.name = data["name"]
        sequence.description = data.get("description", None)
        for a in data["annotations"]:
            sequence.annotate(
                a["start"], a["end"] - 1, a["name"], a["type"], a["color"]
            )
        return sequence

    def _collect_cut_sites(self, enzyme_site, cut1=None, cut2=None):
        if hasattr(enzyme_site, "charac"):
            cut1 = enzyme_site.charac[0]
            cut2 = enzyme_site.charac[1]
            enzyme_site = enzyme_site.charac[4]

        if isinstance(enzyme_site, str):
            enzyme_site = Sequence(enzyme_site)
        cut_sites = []
        for match in self.find_iter(enzyme_site):
            cut_sites.append(match.span[0] + cut1)
            cut_sites.append(match.span[1] + cut2 + 1)
        return cut_sites

    def digest(self, enzymes, as_names=False):
        """Supply either a Bio.RestrictionSite or a tuple of (seq, cut1, cut2)

        e.g. ('GTTTAAAC', 4, -4)

        :param enzymes: either a Bio.RestrictionSite or a tuple of (seq, cut1, cut2)
        :type enzymes: list (of tuple|Bio.RestrictionSite)
        :return: list of sequences
        :rtype: list
        """
        cut_sites = []

        if not isinstance(enzymes, list):
            enzymes = [enzymes]
        if as_names:
            enzymes = [getattr(Restriction, name) for name in enzymes]
        for enzyme in enzymes:
            if isinstance(enzyme, tuple):
                cut_sites += self._collect_cut_sites(*enzyme)
            else:
                cut_sites += self._collect_cut_sites(enzyme)
        return self.cut(cut_sites)

    def __repr__(self) -> str:
        max_width = 30
        replace = "..."
        display = int((max_width - len(replace)) / 2.0)
        s = str(self)
        if len(s) > display * 2:
            # diff = display*2 - len(s)
            s = s[:display] + "..." + s[-display:]
        return "Sequence('{}')".format(s)


Sequence.IO = IOInterface(Sequence)
Sequence.Align = AlignInterface(Sequence)
