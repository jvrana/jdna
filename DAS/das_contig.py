'''
Project: jdna
File: contig
Author: Justin
Date: 2/6/17

Description: 

'''

from collections import defaultdict
import itertools
import json
from copy import deepcopy, copy
import base64
import tempfile
import os
import zipfile
from das_seqio import *
from das_utilities import *


class ContigError(Exception):
    def __init__(self, msg):
        self.msg = msg


class RegionError(Exception):
    def __init__(self, msg):
        self.msg = msg


class ContigContainerMeta(object):
    def __init__(self, source=None, blastver=None, query='untitled', query_length=0, query_circular=False,
                 query_seq=None, contig_seqs=None, **kwargs):
        self.source = source
        self.blastver = blastver
        self.query = query
        self.query_length = query_length
        self.query_circular = query_circular
        self.query_seq = query_seq
        self.contig_seqs = contig_seqs
        self.__dict__.update(kwargs)


class Region(object):
    """
    Classifies a region of a sequence. A region is defined by the inclusive "start" and "end"
    positions in context of an arbitrary sequence defined by the start_index and length.

    Regions can be circular or linear. For circular regions, negative indicies and indicies greater
    than the length are allowable and will be converted to appropriate indices.

    Direction of the region can be FORWARD or REVERSE or BOTH. For reversed directions, start
    and end positions should be flipped.

    Alternative start_index can be used (DEFAULT: 1) to handle sequences that start at 0 or 1.

    A new Region can be created either by defining the start and end positions or defining the
    start position and region_span by Region.create(length, circular)

    E.g. Linear Region
        length: 9
        start_index: 1
        start: 2
        end: 5
        region_span = 5-2+1 = 4
        Context:  |-------|
        C_Index:  1.......9
        Region:    2..4

    E.g. Circular Region
        length: 9
        start_index: 1
        start: 8
        end: 2
        region_span = 4
        Context:  |-------|
        C_Index:  1.......9
        Region:   .2     8.
    """
    START_INDEX = 1
    FORWARD = 1
    REVERSE = -1
    BOTH = 2

    def __init__(self, start, end, length, circular, direction=FORWARD, name=None, start_index=START_INDEX):
        """

        :param start:
        :param end:
        :param length:
        :param circular:
        :param start_index:
        """

        self.__length = length  # not allowed to reset length
        self.__circular = circular  # not allowed to reset length
        self.__bounds_start = start_index  # not allowed to reset length
        self.__start = self.translate_pos(start)
        self.__end = self.translate_pos(end)
        if direction == Region.REVERSE:
            self.__start, self.__end = self.__end, self.__start
        self.__direction = direction  # 1 or -1
        if self.direction not in [Region.FORWARD, Region.REVERSE, Region.BOTH]:
            raise RegionError("Direction {} not understood. Direction must be Region.FORWARD = {}, Region.REVERSE = {},\
             or Region.BOTH = {}".format(self.direction, Region.FORWARD, Region.BACKWARD, Region.BOTH))
        self.name = name
        if not self.circular and self.start > self.end:
            raise RegionError("START cannot be greater than END for linear regions.")

    @staticmethod
    def create(length, circular, direction=FORWARD, name=None, start_index=START_INDEX):
        r = Region(start_index, start_index + 1, length, circular, direction=direction, name=name,
                   start_index=start_index)
        return r

    @property
    def length(self):
        return self.__length

    @property
    def circular(self):
        return self.__circular

    @property
    def bounds_start(self):
        return self.__bounds_start

    @property
    def bounds_end(self):
        return self.length + self.bounds_start - 1

    @property
    def start(self):
        """
        Gets the start position of the region. Internally
        reverses the start and end positions if direction is
        reversed.
        :return:
        """
        return self.__start

    @start.setter
    def start(self, x):
        """
        Sets the start position of the region. Internally
        reverses the start and end positions if direction is
        reversed.
        :return:
        """
        self.__start = self.translate_pos(x)

    @property
    def end(self):
        """
        Gets the end position of the region. Internally
        reverses the start and end positions if direction is
        reversed.
        :return:
        """
        return self.__end

    @end.setter
    def end(self, x):
        """
        Sets the end position of the region. Internally
        reverses the start and end positions if direction is
        reversed.
        :return:
        """
        self.__end = self.translate_pos(x)

    @property
    def direction(self):
        return self.__direction

    def set_region_by_start_and_span(self, x, span):
        if span <= 0:
            raise RegionError("Cannot have a span of less than or equal to zero.")
        self.start = x
        self.end = x + span - 1
        return self

    def set_region_by_delta_start_and_span(self, deltax, span):
        self.set_region_by_start_and_span(deltax + self.bounds_start, span)
        return self

    def translate_pos(self, pos):
        if self.circular:
            cleared = False
            while not cleared:
                cleared = True
                if pos > self.bounds_end:
                    pos = pos - self.length
                    cleared = False
                if pos < self.bounds_start:
                    pos = pos + self.length
                    cleared = False
        else:
            if not self.within_bounds(pos, inclusive=True):
                raise RegionError(
                    "Position {} outside of bounds for linear region [{} {}].".format(pos, self.bounds_start,
                                                                                      self.bounds_end))
        return pos

    @property
    def region_span(self):
        if self.end >= self.start:
            return self.end - self.start + 1
        else:
            if self.circular:
                left = self.bounds_end - self.start + 1
                right = self.end - self.bounds_start + 1
                return left + right
            else:
                raise RegionError("START is greater than END for linear region.")

    def within_region(self, pos, inclusive=True):
        s = self.start
        e = self.end
        if not (self.within_bounds(s, inclusive=True) and self.within_bounds(e, inclusive=True)):
            return False
        if s > e:
            if self.circular:
                if inclusive:
                    return s <= pos <= self.bounds_end or \
                           self.bounds_start <= pos <= e
                else:
                    return s < pos <= self.bounds_end or \
                           self.bounds_start <= pos < e
            else:
                raise RegionError("Cannot have START greater than END for non-circular regions.")
        if inclusive:
            return s <= pos <= e
        else:
            return s < pos < e

    def within_bounds(self, pos, inclusive=True):
        if inclusive:
            return self.bounds_start <= pos <= self.bounds_end
        else:
            return self.bounds_start < pos < self.bounds_end

    def sub_region(self, s, e):
        if self.within_region(s, inclusive=True) and self.within_region(e, inclusive=True):
            r = Region(s, e, self.length, self.circular, direction=self.direction, name=self.name,
                       start_index=self.bounds_start)
            return r
        else:
            raise RegionError(
                "Sub region bounds [{}-{}] outside of Region bounds [{}-{}]".format(s, e, self.bounds_start,
                                                                                    self.bounds_end))

    def same_context(self, other):
        return self.circular == other.circular and \
               self.bounds_start == other.bounds_start and \
               self.bounds_end == other.bounds_end

    def copy(self):
        s, e = self.start, self.end
        if self.direction is Region.REVERSE:
            s, e = self.end, self.start
        return Region(s, e, self.length, self.circular,
                      direction=self.direction, name=self.name, start_index=self.bounds_start)

    def get_overlap(self, other):

        if self.end_overlaps_with(other):
            r = self.copy()
            r.start = other.start
            r.end = self.end

            return r
        else:
            return None

    def end_overlaps_with(self, other):
        """
         Whether this region overlaps the next region it this regions end
             True
                 self   |------|
                 other      |-------|

             False
                 self         |------|
                 other  |-------|

             False
                 self   |------|
                 other    |----|
         :param other:
         :return:
         """

        return self.same_context(other) \
               and self.within_region(other.start, inclusive=True) \
               and not self.within_region(other, inclusive=True)


    def get_gap(self, other):
        if self.consecutive_with(other):
            return None
        if self.no_overlap(other):
            r = self.copy()
            r.start = self.end + 1
            r.end = other.start - 1
            return r
        else:
            return None

    def get_gap_degree(self, other):
        overlap = self.get_overlap(other)
        gap = self.get_gap(other)
        cons = self.consecutive_with(other)
        if cons:
            return 0
        if overlap is not None:
            return -overlap.region_span
        if gap is not None:
            return gap.region_span

    def no_overlap(self, other):
        return self.same_context(other) \
                and not self.within_region(other.start, inclusive=True) \
                and not other.within_region(self.end, inclusive=True)


    def consecutive_with(self, other):
        after_self_pos = None
        before_other_pos = None
        try:
            before_other_pos = other.translate_pos(other.start - 1)
        except RegionError as e:
            return False
        try:
            after_self_pos = self.translate_pos(self.end + 1)
        except RegionError as e:
            return False
        fusable = self.same_context(other) and \
                  other.start == after_self_pos and \
                  self.end == before_other_pos
        return fusable


    def fuse(self, other):
        if self.consecutive_with(other):
            self.end += other.region_span
            return self
        else:
            raise RegionError(
                "Cannot fuse regions [{}-{}] with [{}-{}].".format(self.start, self.end, other.start, other.end))


class Contig(object):
    contig_id = 0
    NEW_PRIMER = "new_primer"
    DIRECT_END = "direct"

    START_INDEX = 1  # Convention is carried over from BLAST results, BE CAREFUL!

    DEFAULT_END = NEW_PRIMER

    TYPE_BLAST = 'blast'
    # TYPE_DIRECT = 'direct_synthesis'
    TYPE_PRIMER = 'primer'
    TYPE_PCR = 'product'
    sequence_options = ['coral', 'seqrecord']
    optional_args = ['circular', 'filename'] + sequence_options

    @staticmethod
    def create_default_contig():
        return Contig(
            q_start=Contig.START_INDEX,
            q_end=Contig.START_INDEX,
            s_end=Contig.START_INDEX,
            s_start=Contig.START_INDEX,
            identical=0,
            query_acc='',
            subject_acc='',
            query_seq='',
            contig_type='',
            start_label=Contig.DEFAULT_END,
            end_label=Contig.DEFAULT_END,
            gaps=0,
            bit_score=0,
            quality=0,
            subject_strand='',
            evalue=0,
            gap_opens=0,
            filename='',
            subject_seq='',
            score=0,
            subject_length=0,
            query_length=0,
            alignment_length=0,
            circular=False
        )

    def __init__(self, **kwargs):
        subject_direction = Region.FORWARD
        if kwargs['subject_strand'] == 'minus':
            subject_direction = Region.REVERSE
        # Subject region information
        self.subject = Region(
            kwargs['s_start'],
            kwargs['s_end'],
            kwargs['subject_length'],
            kwargs['subject_circular'],
            direction=subject_direction,
            name=kwargs['subject_acc'],
            start_index=Contig.START_INDEX
        )
        self.subject.strand = kwargs['subject_strand']
        self.subject.seq = kwargs['subject_seq']

        # Query region information
        self.query = Region(
            kwargs['q_start'],
            kwargs['q_end'],
            kwargs['query_length'],
            kwargs['query_circular'],
            name=kwargs['query_acc'],
            direction=Region.FORWARD,
            start_index=Contig.START_INDEX
        )
        self.query.seq = kwargs['query_seq']

        # Alignment Information
        self.score = kwargs['score']
        self.evalue = kwargs['evalue']
        self.bit_score = kwargs['bit_score']
        self.alignment_length = kwargs['alignment_length']  # for comparison
        self.identical = kwargs['identical']
        self.gap_opens = kwargs['gap_opens']
        self.gaps = kwargs['gaps']

        self.contig_type = kwargs['contig_type']
        self.assign_id()
        if 'end_label' not in kwargs:
            self.end_label = Contig.DEFAULT_END
        else:
            self.end_label = kwargs['end_label']

        if 'start_label' not in kwargs:
            self.start_label = Contig.DEFAULT_END
        else:
            self.start_label = kwargs['start_label']

        self.quality = 1.0  # used for future machine learning / AI

        for opt in Contig.optional_args:
            try:
                self.__dict__[opt] = kwargs[opt]
            except KeyError:
                self.__dict__[opt] = None

        # TODO: this implies a direct synthesis, change this to add separate
        if not self.circular and self.is_perfect_subject():
            self.start_label = Contig.DIRECT_END
            self.end_label = Contig.DIRECT_END
            # for k in kwargs:
            #     if k not in self.__dict__:
            #         raise ValueError("Key {} not found in {} class definition".format(k, Contig.__class__.__name__))

    def assign_id(self):
        Contig.contig_id += 1
        self.contig_id = Contig.contig_id

    def deepcopy(self):
        c = deepcopy(self)
        c.assign_id()
        return c

    def is_direct(self):
        return self.start_label == Contig.DIRECT_END and self.end_label == Contig.DIRECT_END

    def assign_id(self):
        Contig.contig_id += 1
        self.contig_id = Contig.contig_id

    def copy(self):
        c = deepcopy(self)
        c.assign_id()
        return c

    def divide_contig(self, startpoints, endpoints, include_contig=True, contig_type=None, circular=None):
        """Divides contigs by start and end points. Contigs inherit end_labels. If end labels are None,
        default value of Contig.NEW_PRIMER is used. Ignores positions outside of bounds.

        e.g.
        OLD                |--------------------------------|
        Points                (x)>      <   >
                           |-----|
                              (x)|------|
                                        |---|
                                            |----------------|
                                        |--------------------|
                              (x)|---------------------------|
                           |---------------------------------|     <<< if include_contig == True
                              (x)|----------|
                           |------------|
        :param startpoints: list of (pos, label) tuples
        :param endpoints: list of (pos, label) tuples
        :param include_contig: whether to include the contig endpoints
        :param contig_type: contig_label, used to determine the type of contig this is for the cost analysis
        :param circular: used to determine the type of contig this is for the cost analysis
        :return:
        """

        positions = list(itertools.product(startpoints, endpoints))
        contigs = []

        for start, end in positions:
            s, start_label = start
            e, end_label = end
            if not include_contig:
                if s == self.query.start and e == self.query.end:
                    continue
            try:
                new_contig = self.break_contig(s, e, start_label=start_label, end_label=end_label)
                if contig_type is not None:
                    new_contig.contig_type = contig_type
                if circular is not None:
                    new_contig.circular = circular
                new_contig.parent_id = self.contig_id
                contigs.append(new_contig)
            except ContigError as e:
                pass
        return contigs

    def equivalent_location(self, other):
        return other.query.start == self.query.start and other.query.end == self.query.end

    def overlaps(self, other, inclusive=True):
        return other.query.within_region(self.query.start, inclusive=inclusive) or \
               other.query.within_region(self.query.end, inclusive=inclusive)

    def is_within(self, other, inclusive=True):
        return other.query.within_region(self.query.start, inclusive=inclusive) and \
               other.query.within_region(self.query.end, inclusive=inclusive)

    # def is_within(self, other, inclusive=True):
    #     return self.query.start >= other.query.start and self.query.end <= other.query.end
    # TODO: Reevaluate subject start and subject end after breaking
    # TODO: Circular break contig; but its much faster to do it the way we are now...
    def break_contig(self, q_start, q_end, start_label=None, end_label=None):
        """
        Breaks a contig at query start and end; copies over information
        from self contig
        :param q_start:
        :param q_end:
        :param start_label:
        :param end_label:
        :return:
        """
        if q_start > q_end:
            raise ContigError("query_start cannot be greater than query_end")
        if not (self.query.within_region(q_start, inclusive=True) and self.query.within_region(q_end, inclusive=True)):
            e = "break points [{}, {}] are outside bounds of contig bounds [{}, {}]".format(q_start, q_end,
                                                                                            self.query.start,
                                                                                            self.query.end)
            raise ContigError('msg')
        new_contig = self.deepcopy()

        # Set query region
        new_contig.query.start = q_start
        new_contig.query.end = q_end

        # Set subject region
        delta_start = q_start - self.query.start
        delta_end = q_end - self.query.end
        new_contig.subject.start += delta_start
        new_contig.subject.end += delta_end

        new_contig.alignment_length = new_contig.query.region_span
        new_contig.parent_id = self.contig_id
        if start_label is not None:
            new_contig.start_label = start_label
        else:
            if new_contig.query.start == self.query.start:
                new_contig.start_label = self.start_label
            else:
                new_contig.start_label = Contig.DEFAULT_END
        if end_label is not None:
            new_contig.end_label = end_label
        else:
            if new_contig.query.end == self.query.end:
                new_contig.end_label = self.end_label
            else:
                new_contig.end_label = Contig.DEFAULT_END
        return new_contig

    def is_perfect_subject(self):
        return self.alignment_length == self.subject.length and self.is_perfect_alignment()

    def is_perfect_alignment(self):
        return self.alignment_length == self.identical and \
               self.gaps == 0 and self.gap_opens == 0

    def pcr_products_of_contig(self, primers):
        """

        :param primers:
        :return:
        """
        primers = self.get_primers_within_bounds(primers)

        rev_primer_pos = []
        fwd_primer_pos = []
        new_fwd_primer_pos = [(self.query.start, None)]
        new_rev_primer_pos = [(self.query.end, None)]
        for primer in primers:
            direction = primer.subject.strand
            if direction == 'plus':
                if self.query.within_region(primer.query.start):
                    fwd_primer_pos.append((primer.query.start, primer.contig_id))
                    new_rev_primer_pos.append((primer.query.end, None))
            if direction == 'minus':
                if self.query.within_region(primer.query.end):
                    rev_primer_pos.append((primer.query.end, primer.contig_id))
                    new_fwd_primer_pos.append((primer.query.start, None))
        contigs = []
        contigs += self.divide_contig(fwd_primer_pos, rev_primer_pos, include_contig=False,
                                      contig_type=Contig.TYPE_PCR, circular=False)
        contigs += self.divide_contig(fwd_primer_pos, new_rev_primer_pos, include_contig=False,
                                      contig_type=Contig.TYPE_PCR, circular=False)
        contigs += self.divide_contig(new_fwd_primer_pos, rev_primer_pos, include_contig=False,
                                      contig_type=Contig.TYPE_PCR, circular=False)
        contigs += self.divide_contig(new_fwd_primer_pos, new_rev_primer_pos, include_contig=False,
                                      contig_type=Contig.TYPE_PCR, circular=False)
        return contigs

    def get_primers_within_bounds(self, primers, minimum_primer_anneal=15):

        def primer_within_bounds(primer):
            if primer.subject.strand == 'plus':
                return self.query.start + minimum_primer_anneal < primer.query.end < self.query.end
            elif primer.subject.strand == 'minus':
                return self.query.start < primer.query.start < self.query.end - minimum_primer_anneal
            return False

        return filter(lambda x: primer_within_bounds(x), primers)

    def ends_equivalent(self, other):
        return self.start_label == other.start_label and \
               self.end_label == other.end_label

    def json(self):
        j = self.__dict__
        for x in Contig.sequence_options:
            del j[x]
        j['query'] = self.query.__dict__
        j['subject'] = self.subject.__dict__
        return j


class PerfectPrimerContig(Contig):
    """
    A perfect contig with no gaps (i.e. a perfect primer)
    """

    def __init__(self, query_acc, q_start, q_end, subject_acc, s_start, s_end, strand, filename):
        super(PerfectPrimerContig, self).__init__(
            query_acc=query_acc,
            subject_acc=subject_acc,
            q_start=q_start,
            q_end=q_end,
            s_start=s_start,
            s_end=s_end,
            identical=s_end - s_start + 1,
            contig_type="perfect_contig",
            start_label=None,
            end_label=None,
            gaps=0,
            bit_score=0.0,
            quality=1.0,
            subject_strand=strand,
            evalue=0.0,
            gap_opens=0,
            filename=filename,
            subject_seq=''
        )


class ContigContainer(object):
    def __init__(self, meta=None, contigs=None):
        if meta is None:
            meta = {}
        self.add_metadata(**meta)
        if contigs is None:
            contigs = []
        self.contigs = contigs
        self.make_dictionary()
        # self.sort_contigs()

    def make_dictionary(self):
        self.contig_dictionary = {x.contig_id: x for x in self.contigs}

    def get_contig(self, contig_id):
        return self.contig_dictionary[contig_id]

    def add_contig(self, **kwargs):
        new_contig = Contig(**kwargs)
        self.contigs.append(new_contig)
        self.make_dictionary()
        return new_contig

    def add_metadata(self, **kwargs):
        if kwargs is None:
            kwargs = {}
        self.meta = ContigContainerMeta(**kwargs)

    # TODO: add alternative templates field to contig
    def remove_redundant_contigs(self, remove_equivalent=True, remove_within=False, no_removal_if_different_ends=True):
        """
        Selects contigs based on criteria
        :param remove_equivalent:
        :param remove_within:
        :param no_removal_if_different_ends:
        :return:
        """

        contigs_for_removal = []
        for c1 in self.contigs:
            for c2 in self.contigs:
                if c1 == c2:
                    continue
                if c1 in contigs_for_removal or c2 in contigs_for_removal:
                    continue
                #
                if no_removal_if_different_ends:
                    if not c1.ends_equivalent(c2):
                        continue
                if remove_equivalent and c1.equivalent_location(c2):
                    contigs_for_removal.append(c2)
                if remove_within:
                    if c1.is_within(c2, inclusive=True):
                        contigs_for_removal.append(c1)
                    elif c2.is_within(c1, inclusive=True):
                        contigs_for_removal.append(c2)
        contigs_for_removal = list(set(contigs_for_removal))
        for c in contigs_for_removal:
            self.contigs.remove(c)
        self.make_dictionary()
        # TODO: prefer to remove contigs with new primers over contigs with existing primers

    def fuse_circular_fragments(self):
        """
        Fuses contigs that are circular and would
        span the origin if the index of the blast subject
        had been different
        :return: None
        """
        ga = defaultdict(list)

        def fuse_condition(l, r):
            return l.subject.consecutive_with(r.subject) and \
                   l.query.consecutive_with(r.query) and \
                   l.subject.name == r.subject.name and \
                   l in self.contigs and \
                   r in self.contigs

        def fuse(l, r):
            l.query.fuse(r.query)
            l.subject.fuse(r.subject)

            keys_to_sum = 'alignment_length, gap_opens, gaps, identical, score'.split(', ')
            for k in keys_to_sum:
                l.__dict__[k] += r.__dict__[k]

        pairs = itertools.permutations(self.contigs, 2)
        for l, r in pairs:
            if fuse_condition(l, r):
                fuse(l, r)
                self.contigs.remove(r)

    def dump(self, out):
        """
        Dumps this container to a json
        :param out:
        :return:
        """
        j = deepcopy(self.__dict__)
        del j['contig_dictionary']
        j['meta'] = j['meta'].__dict__
        j['contigs'] = [x.json() for x in j['contigs']]
        del j['meta']['query_seq']
        del j['meta']['contig_seqs']

        with open(out, 'w') as output:
            json.dump(j, output)

    def sort_contigs(self):
        """
        Sorts contigs according to their query start
        :return:
        """
        self.contigs = sorted(self.contigs, key=lambda x: x.query.start)

    def filter_perfect_subjects(self):
        """
        Selects only contigs with perfect subjects, or with
        subjects that have 100% of their identity contained
        in the contig (no partial subjects)
        :return: None
        """
        filtered_contigs = []
        for contig in self.contigs:
            if contig.is_perfect_subject():
                filtered_contigs.append(contig)
        self.contigs = filtered_contigs

    def filter_perfect(self):
        """
        Selects only contigs that have perfect alignment scores and are perfect subjects
        i.e. no gaps, no mismatches, 100% of subjects
        :return:
        """
        filtered_contigs = []
        for contig in self.contigs:
            if contig.is_perfect_alignment():
                filtered_contigs.append(contig)
            else:
                pass
        self.contigs = filtered_contigs
        # TODO: handle ambiquoous NNNN dna in blast search by eliminating gap_opens, gaps if they are N's

    def expand_contigs(self, primers):
        """
        Creates PCR products for all contigs from a container of primers
        e.g. with "P>" = Primer
        Primer        P>   <P
        OLD  |--------------------|
        NEW           P|----------|   (new reverse primer)
        NEW  |--------|               (two new primers)
        NEW           P|----|P        (uses forward and reverse primers)
        NEW  |--------------|P        (new forward primer)
        :param primers: primer container
        :return:
        """
        all_alignments = []
        for contig in self.contigs:
            new_alignments = contig.pcr_products_of_contig(primers)
            all_alignments += new_alignments
        self.contigs += all_alignments

    # def break_long_contigs(self):
    #     """
    #     Breaks up contigs whose length is greater than the allowable subject length
    #     :return:
    #     """
    #     new_contigs = []
    #     for c in self.contigs:
    #         if c.get_length


    def break_contigs_at_endpoints(self, contig_type=Contig.TYPE_PCR):
        """
        Breaks long contigs at endpoints of overlapping contigs
        e.g. with "X"=break point & "|" = endpoint
            OLD   |---------------X-------|
            OLD                   |-----------------|
            NEW   |---------------|
            NEW                   |-------|
        :param contig_type: type label of contig to assign newly created contigs
        :return: None
        """
        new_contigs = []
        for c1 in self.contigs:
            for c2 in self.contigs:
                if c1.query.within_region(c2.query.end, inclusive=False):  # if second contig overlaps with first contig
                    # new_contigs.append(c1.break_contig(c1.query.start, c2.query.end, start_label=None, end_label=None))
                    new_contigs.append(c1.break_contig(c2.query.end, c1.query.end, start_label=None, end_label=None))
                if c1.query.within_region(c2.query.start, inclusive=False):
                    new_contigs.append(
                        c1.break_contig(c1.query.start, c2.query.start, start_label=None, end_label=None))
                    # new_contigs.append(c1.break_contig(c2.query.start, c1.query.end, start_label=None, end_label=None))
        for n in new_contigs:
            n.contig_type = contig_type + ' (broken long contig)'

        for c1 in self.contigs:
            end = c1.query.start + self.meta.query_length
            if end < c1.query.length:
                for c2 in self.contigs:
                    if c2.query.start < end < c2.query.end:
                        n = c2.break_contig(c2.query.start, end)
                        n.contig_type = contig_type + ' (broken for circularization)'
                        new_contigs.append(n)

        self.contigs += new_contigs
