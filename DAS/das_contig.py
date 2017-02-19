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
import xmlrpclib
from das_seqio import *
from das_utilities import *

class ContigError(Exception):
    def __init__(self, msg):
        self.msg = msg


class ContigContainerMeta(object):
    def __init__(self, source=None, blastver=None, query='untitled', query_length=0, query_circular=False,
                 query_seq=None, contig_seqs=None, **kwargs):
        self.source = source
        self.blastver = blastver
        self.query = query
        self.query.length = query_length
        self.query_circular = query_circular
        self.query_seq = query_seq
        self.contig_seqs = contig_seqs
        self.__dict__.update(kwargs)

class Region(object):
    """
    Classifies a region of DNA; can be circular or linear
    """
    START_INDEX = 1



    def __init__(self, start, end, length, circular, name=None, start_index=START_INDEX):
        """

        :param start:
        :param end:
        :param length:
        :param circular:
        :param start_index:
        """
        self.start = self.translate_pos(start)
        self.end = self.translate_pos(end)
        self.length = length
        self.circular = circular
        self._start_index = start_index
        self.name = name
        self.verify()


    def verify(self):
        assert self._start_index <= self.start <= self.length + self._start_index
        assert self._start_index <= self.end <= self.length + self._start_index

    def translate_pos(self, pos):
        if self.circular:
            cleared = False
            while not cleared:
                cleared = True
                if pos >= self.length + self._start_index:
                    pos = pos - self.length
                    cleared = False
                if pos < self._start_index:
                    pos = pos + self.length
                    cleared = False
        else:
            assert self.within(pos)
        return pos

    def within(self, pos, inclusive=True):
        if inclusive:
            return self.get_start() <= pos <= self.get_end()
        else:
            return self.get_start() < pos < self.get_end()

    def get_span(self):
        return self.end - self.start + 1

    def get_start(self):
        return self._start_index

    def get_end(self):
        return self.length + self.get_start() - 1



class Contig(object):
    contig_id = 0
    NEW_PRIMER = "new_primer"
    DIRECT_END = "direct"

    START_INDEX = 1 # Convention is carried over from BLAST results, BE CAREFUL!

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
        super(Contig, self).__init__(**kwargs)

        # Subject region information
        self.subject = Region(
            kwargs['s_start'],
            kwargs['s_end'],
            kwargs['subject_length'],
            kwargs['circular'],
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
            kwargs['circular'],
            name=kwargs['query_acc'],
            start_index=Contig.START_INDEX
        )
        self.query.seq = kwargs['query_seq']

        # Alignment Information
        self.score = kwargs['score']
        self.evalue = kwargs['evalue']
        self.bit_score = kwargs['bit_score']
        self.alignment_length = kwargs['alignment_length'] # for comparison
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
        if not self.circular:
            self.start_label = Contig.DIRECT_END
            self.end_label = Contig.DIRECT_END
        for k in kwargs:
            if k not in self.__dict__:
                raise ValueError("Key {} not found in {} class definition".format(k, Contig.__class__.__name__))

    def json(self):
        return self.__dict__

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

    def q_pos_within(self, pos, inclusive=True):
        if inclusive:
            return self.query.start <= pos <= self.query.end
        else:
            return self.query.start < pos < self.query.end

    def overlaps(self, other, inclusive=True):
        return other.q_pos_within(self.query.start, inclusive=inclusive) or \
               other.q_pos_within(self.query.end, inclusive=inclusive)

    def is_within(self, other, inclusive=True):
        return other.q_pos_within(self.query.start, inclusive=inclusive) and \
               other.q_pos_within(self.query.end, inclusive=inclusive)

    # def is_within(self, other, inclusive=True):
    #     return self.query.start >= other.query.start and self.query.end <= other.query.end

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
        if not (self.q_pos_within(q_start, inclusive=True) and self.q_pos_within(q_end, inclusive=True)):
            e = "break points [{}, {}] are outside bounds of contig bounds [{}, {}]".format(q_start, q_end, self.query.start, self.query.end)
            raise ContigError('msg')
        new_contig = self.deepcopy()
        new_contig.query.start = q_start
        new_contig.query.end = q_end



        new_contig.subject.start = convert_circular_position(self.subject.start + (q_start - self.query.start), self.subject.length, 1)
        new_contig.subject.end = convert_circular_position(self.subject.end - (self.query.end - q_end), self.subject.length, 1)



        new_contig.alignment_length = q_end - q_start
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
            direction = primer.subject_strand
            if direction == 'plus':
                if self.q_pos_within(primer.query.start):
                    fwd_primer_pos.append((primer.query.start, primer.contig_id))
                    new_rev_primer_pos.append((primer.query.end, None))
            if direction == 'minus':
                if self.q_pos_within(primer.query.end):
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
            if primer.subject_strand == 'plus':
                return self.query.start + minimum_primer_anneal < primer.query.end < self.query.end
            elif primer.subject_strand == 'minus':
                return self.query.start < primer.query.start < self.query.end - minimum_primer_anneal
            return False

        return filter(lambda x: primer_within_bounds(x), primers)


    def ends_equivalent(self, other):
        return self.start_label == other.start_label and \
                self.end_label == other.end_label

    def json(self):
        j = super(Contig, self).json()
        for x in Contig.sequence_options:
            del j[x]
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
            identical=s_end-s_start+1,
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
            return l.subject.name == r.subject.name and \
                   l.query.end + 1 == r.query.start and \
                   l.subject.length == l.subject.end and \
                   l.circular and r.circular and \
                    l in self.contigs and \
                    r in self.contigs

        def fuse(l, r):
            l.query.end = r.query.end
            l.subject.end = r.subject.end

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
                if c1.q_pos_within(c2.query.end, inclusive=False): # if second contig overlaps with first contig
                    # new_contigs.append(c1.break_contig(c1.query.start, c2.query.end, start_label=None, end_label=None))
                    new_contigs.append(c1.break_contig(c2.query.end, c1.query.end, start_label=None, end_label=None))
                if c1.q_pos_within(c2.query.start, inclusive=False):
                    new_contigs.append(c1.break_contig(c1.query.start, c2.query.start, start_label=None, end_label=None))
                    # new_contigs.append(c1.break_contig(c2.query.start, c1.query.end, start_label=None, end_label=None))
        for n in new_contigs:
            n.contig_type = contig_type + ' (broken long contig)'

        for c1 in self.contigs:
            end = c1.query.start + self.meta.query.length
            if end < c1.query.length:
                for c2 in self.contigs:
                    if c2.query.start < end < c2.query.end:
                        n = c2.break_contig(c2.query.start, end)
                        n.contig_type = contig_type + ' (broken for circularization)'
                        new_contigs.append(n)

        self.contigs += new_contigs
