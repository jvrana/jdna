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
from copy import deepcopy

class ContigError(Exception):

    def __init__(self, message):
        self.message = message

class ContigContainerMeta(object):

    def __init__(self, source=None, blastver=None, query='untitled', query_length=0, query_circular=False, **kwargs):
        self.source = source
        self.query = query
        self.query_length = query_length
        self.query_circular = query_circular
        self.__dict__ = kwargs

class ContigContainer(object):

    def __init__(self, meta=None, contigs=[]):
        if meta is None:
            meta = {}
        self.add_metadata(**meta)
        self.contigs = contigs

    def make_dictionary(self):
        return {x.contig_id: x for x in self.contigs}

    def get_contig(self, id):
        return self.make_dictionary()[id]

    def add_contig(self, **kwargs):
        new_contig = Contig(**kwargs)
        self.contigs.append(new_contig)
        return new_contig

    def add_metadata(self, **kwargs):
        if kwargs is None:
            kwargs = {}
        self.meta = ContigContainerMeta(**kwargs)

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
            keys_to_sum = 'alignment_length, identical, gap_opens, gaps, identical, score'.split(', ')
            for k in keys_to_sum:
                l.__dict__[k] += r.__dict__[k]

        pairs = itertools.permutations(self.contigs, 2)
        for l, r in pairs:
            if fuse_condition(l, r):
                fuse(l, r)
                self.contigs.remove(r)

    def dump(self, out):
        j = deepcopy(self.__dict__)
        j['meta'] = j['meta'].__dict__
        j['contigs'] = [x.__dict__ for x in j['contigs']]
        with open(out, 'w') as output:
            json.dump(j, output)


class Contig(object):
    contig_id = 0

    def __init__(self, **kwargs):
        self.query_acc = kwargs['query_acc']
        self.subject_acc = kwargs['subject_acc']
        self.score = kwargs['score']
        self.evalue = kwargs['evalue']
        self.bit_score = kwargs['bit_score']
        self.alignment_length = kwargs['alignment_length']
        self.identical = kwargs['identical']
        self.gap_opens = kwargs['gap_opens']
        self.gaps = kwargs['gaps']
        self.query_length = kwargs['query_length']
        self.q_start = kwargs['q_start']
        self.q_end = kwargs['q_end']
        self.subject_length = kwargs['subject_length']
        self.s_start = kwargs['s_start']
        self.s_end = kwargs['s_end']
        self.subject_strand = kwargs['subject_strand']
        self.query_seq = kwargs['query_seq']
        self.subject_seq = kwargs['subject_seq']
        self.assemblies = []
        self.contig_type = kwargs['contig_type']
        self.assign_id()
        try:
            self.circular = kwargs['circular']
        except KeyError:
            self.circular = False
        try:
            self.filename = kwargs['filename']
        except KeyError:
            self.filename = None
        for k in kwargs:
            if k not in self.__dict__:
                raise ValueError("Key {} not found in {} class definition".format(k, Contig.__class__.__name__))

    def assign_id(self):
        Contig.contig_id += 1
        self.contig_id = Contig.contig_id

    def copy(self):
        c = deepcopy(self)
        c.assign_id()
        return c

    def divide_contig(self, startpoints, endpoints, includecontig=True, contig_type=None):
        if includecontig:
            startpoints += [self.q_start]
            endpoints += [self.q_end]
        positions = list(itertools.product(startpoints, endpoints))
        contigs = []
        for s, e in positions:
            try:
                new_contig = self.break_contig(s, e)
                contigs.append(new_contig)
                if contig_type is not None:
                    new_contig.contig_type = contig_type
            except ContigError:
                pass
        return contigs

    def break_contig(self, q_start, q_end):
        '''
        Breaks a contig at query start and end; copies over information
        from self contig
        :param q_start:
        :param q_end:
        :return:
        '''
        if q_start < self.q_start:
            raise ContigError("query_start {} cannot be less than contig start {}".format(q_start, self.q_start))
        if q_end > self.q_end:
            raise ContigError("query_end {} cannot be less than contig end {}".format(q_start, self.q_start))
        if q_start > q_end:
            raise ContigError("query_start cannot be greater than query_end")
        new_contig = self.copy()
        new_contig.q_start = q_start
        new_contig.q_end = q_end
        new_contig.s_start = self.s_start + (q_start - self.q_start)
        new_contig.s_end = self.s_end - (self.q_end - q_end)
        new_contig.alignment_length = q_end - q_start
        return new_contig

