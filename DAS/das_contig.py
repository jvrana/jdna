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

class ContigError(Exception):
    def __init__(self, message):
        self.message = message


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


class QueryRegion(object):
    contig_id = 0

    def __init__(self, **kwargs):
        self.query_acc = kwargs['query_acc']
        self.q_start = kwargs['q_start']
        self.q_end = kwargs['q_end']
        self.query_length = kwargs['query_length']
        self.assign_id()

    def get_span(self):
        return self.q_end - self.q_start

    def json(self):
        return self.__dict__

    def assign_id(self):
        Contig.contig_id += 1
        self.contig_id = Contig.contig_id

    def deepcopy(self):
        c = deepcopy(self)
        c.assign_id()
        return c


class Contig(QueryRegion):
    NEW_PRIMER = "new_primer"
    DIRECT_END = "direct"

    TYPE_BLAST = 'blast'
    # TYPE_DIRECT = 'direct_synthesis'
    TYPE_PRIMER = 'primer'
    TYPE_PCR = 'product'
    sequence_options = ['coral', 'seqrecord']
    optional_args = ['circular', 'filename'] + sequence_options

    def __init__(self, **kwargs):
        super(Contig, self).__init__(**kwargs)
        self.subject_acc = kwargs['subject_acc']
        self.score = kwargs['score']
        self.evalue = kwargs['evalue']
        self.bit_score = kwargs['bit_score']
        self.alignment_length = kwargs['alignment_length']
        self.identical = kwargs['identical']
        self.gap_opens = kwargs['gap_opens']
        self.gaps = kwargs['gaps']
        self.subject_length = kwargs['subject_length']
        self.s_start = kwargs['s_start']
        self.s_end = kwargs['s_end']
        self.subject_strand = kwargs['subject_strand']
        self.query_seq = kwargs['query_seq']
        self.subject_seq = kwargs['subject_seq']
        self.contig_type = kwargs['contig_type']
        self.assign_id()
        self.end_label = Contig.NEW_PRIMER  # None, new_primer, or contig_id
        self.start_label = Contig.NEW_PRIMER  # None, new_primer, or contig_id
        self.quality = 1.0  # used for future machine learning / AI

        for opt in Contig.optional_args:
            try:
                self.__dict__[opt] = kwargs[opt]
            except KeyError:
                self.__dict__[opt] = None

        # TODO: this implies a direct synthesis
        if not self.circular:
            self.start_label = Contig.DIRECT_END
            self.end_label = Contig.DIRECT_END
        for k in kwargs:
            if k not in self.__dict__:
                raise ValueError("Key {} not found in {} class definition".format(k, Contig.__class__.__name__))

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
        """Divides contigs by start and end points
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
                if s == self.q_start and e == self.q_end:
                    continue
            try:
                new_contig = self.break_contig(s, e)
                if contig_type is not None:
                    new_contig.contig_type = contig_type
                if circular is not None:
                    new_contig.circular = circular
                if new_contig.q_start == self.q_start:
                    new_contig.start_label = Contig.NEW_PRIMER
                else:
                    new_contig.start_label = start_label
                if new_contig.q_end == self.q_end:
                    new_contig.end_label = Contig.NEW_PRIMER
                else:
                    new_contig.end_label = end_label
                new_contig.parent = self.contig_id
                contigs.append(new_contig)
            except ContigError:
                pass
        return contigs

    def equivalent_location(self, other):
        return other.q_start == self.q_start and other.q_end == self.q_end

    def is_within(self, other):
        return self.q_start >= other.q_start and self.q_end <= other.q_end

    def break_contig(self, q_start, q_end):
        """
        Breaks a contig at query start and end; copies over information
        from self contig
        :param q_start:
        :param q_end:
        :return:
        """
        if q_start < self.q_start:
            raise ContigError("query_start {} cannot be less than contig start {}".format(q_start, self.q_start))
        if q_end > self.q_end:
            raise ContigError("query_end {} cannot be less than contig end {}".format(q_start, self.q_start))
        if q_start > q_end:
            raise ContigError("query_start cannot be greater than query_end")
        new_contig = self.deepcopy()
        new_contig.q_start = q_start
        new_contig.q_end = q_end
        new_contig.s_start = self.s_start + (q_start - self.q_start)
        new_contig.s_end = self.s_end - (self.q_end - q_end)
        new_contig.alignment_length = q_end - q_start
        return new_contig

    def is_perfect(self):
        return self.alignment_length == self.identical and \
               self.gaps == 0 and self.gap_opens == 0

    # TODO: add primer label names to start_label and end_label
    def pcr_products_of_contig(self, primers):
        """

        """
        primers = self.get_primers_within_bounds(primers)

        rev_primer_pos = []
        fwd_primer_pos = []
        new_fwd_primer_pos = [(self.q_start, "new_primer")]
        new_rev_primer_pos = [(self.q_end, "new_primer")]
        for primer in primers:
            direction = primer.subject_strand
            if direction == 'plus':
                fwd_primer_pos.append((primer.q_start, primer.contig_id))
                new_rev_primer_pos.append((primer.q_end, "new_primer"))
            if direction == 'minus':
                rev_primer_pos.append((primer.q_end, primer.contig_id))
                new_fwd_primer_pos.append((primer.q_start, "new_primer"))
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
                return self.q_start + minimum_primer_anneal < primer.q_end < self.q_end
            elif primer.subject_strand == 'minus':
                return self.q_start < primer.q_start < self.q_end - minimum_primer_anneal
            return False

        return filter(lambda x: primer_within_bounds(x), primers)

        # TODO: generate additional pcr products that could be homologous to existing primer
        # TODO: compute alignment_graph for each 'contig'

    def json(self):
        j = super(Contig, self).json()
        for x in Contig.sequence_options:
            del j[x]
        return j


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
    def remove_redundant_contigs(self, include_contigs_contained_within=False, save_linear_contigs=False):
        contig_for_removal = []
        for c1 in self.contigs:
            for c2 in self.contigs:
                if c1 == c2:
                    continue
                if c1 in contig_for_removal or c2 in contig_for_removal:
                    continue
                if not c2.circular and save_linear_contigs:
                    continue
                if c1.equivalent_location(c2):
                    contig_for_removal.append(c2)
                if include_contigs_contained_within and c2.is_within(c1):
                    contig_for_removal.append(c2)
        for c in contig_for_removal:
            self.contigs.remove(c)
        self.make_dictionary()

    def fuse_circular_fragments(self):
        ga = defaultdict(list)

        def fuse_condition(l, r):
            return l.subject_acc == r.subject_acc and \
                   l.q_end + 1 == r.q_start and \
                   l.subject_length == l.s_end and \
                   l.circular and r.circular

        def fuse(l, r):
            l.q_end = r.q_end
            l.s_end = r.s_end
            keys_to_sum = 'alignment_length, gap_opens, gaps, identical, score'.split(', ')
            for k in keys_to_sum:
                l.__dict__[k] += r.__dict__[k]

        pairs = itertools.permutations(self.contigs, 2)
        for l, r in pairs:
            if fuse_condition(l, r):
                fuse(l, r)
                self.contigs.remove(r)

    def dump(self, out):
        j = deepcopy(self.__dict__)
        del j['contig_dictionary']
        j['meta'] = j['meta'].__dict__
        j['contigs'] = [x.json() for x in j['contigs']]
        del j['meta']['query_seq']
        del j['meta']['contig_seqs']
        with open(out, 'w') as output:
            json.dump(j, output)

    def sort_contigs(self):
        self.contigs = sorted(self.contigs, key=lambda x: x.q_start)

    def filter_perfect(self):
        filtered_contigs = []
        for contig in self.contigs:
            if contig.is_perfect():
                filtered_contigs.append(contig)
            else:
                pass

                # TODO: handle ambiquoous NNNN dna in blast search by eliminating gap_opens, gaps if they are N's

    def expand_contigs(self, primers):
        all_alignments = []
        for contig in self.contigs:
            new_alignments = contig.pcr_products_of_contig(primers)
            all_alignments += new_alignments
        self.contigs += all_alignments

