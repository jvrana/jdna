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

    def divide_contig(self, startpoints, endpoints, include_contig=True, contig_type=None):
        if include_contig:
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


class ContigContainer(object):
    def __init__(self, meta=None, contigs=None):
        if meta is None:
            meta = {}
        self.add_metadata(**meta)
        if contigs is None:
            contigs = []
        self.contigs = contigs
        # self.sort_contigs()

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

    def sort_contigs(self):
        self.contigs = sorted(self.contigs, key=lambda x: x.q_start)


class Assembly(ContigContainer):
    def __init__(self, meta=None, contigs=[]):
        super(Assembly, self).__init__(meta=meta, contigs=contigs)
        self.assembly_graph = {}
        self.paths = []

    def get_all_assemblies(self, sort=True, place_holder_size=10):
        self.make_assembly_graph(sort=sort)
        self.assemblies = self.dfs_iter(place_holder_size=place_holder_size)
        return self.assemblies

    def compute_assembly_cost(self, contig_path):
        if len(contig_path) == 0:
            return float("Inf")
        first = contig_path[0]

        pairs = zip(contig_path[:-1], contig_path[1:])

        gaps = [first.q_start]

        for l, r in pairs:
            if first.q_start - r.q_end > first.query_length / 2.0:
                continue
            align_gap = r.q_start - r.q_end
            if align_gap < 0:
                align_gap = 0
            gaps.append(align_gap)

        non_zero_gaps = filter(lambda x: x != 0, gaps)
        cost = sum(non_zero_gaps) + 10 * len(non_zero_gaps) + len(contig_path)

        if cost < 0:
            raise Exception("Cost was calculated to be less than zero.")

        return cost

    def assembly_condition(self, left, right):
        r_5prime_threshold = 0
        r_3prime_threshold = 60
        l_pos = left.q_end
        r_pos = right.q_start
        # return r_pos == l_pos + 1 # if its consecutive
        return r_pos > l_pos - r_3prime_threshold and \
               right.q_end > left.q_end

    def make_assembly_graph(self, sort=True):
        graph = {}
        pairs = list(itertools.permutations(self.contigs, 2))
        print '{} pairs to search'.format(len(pairs))
        for l, r in pairs:

            if l == r:
                continue
            if self.assembly_condition(l, r):
                if l.contig_id not in graph:
                    graph[l.contig_id] = []
                graph[l.contig_id].append(r.contig_id)

        if sort:
            for k in graph.keys():
                a_array = graph[k][:]

                graph[k] = sorted(a_array, key=lambda x: self.compute_assembly_cost(
                    [self.get_contig(k), self.get_contig(x)]))
        self.assembly_graph = graph
        # self.assembly_graph['root'] = [x.contig_id for x in self.contigs]
        return self.assembly_graph

    def dfs_iter(self, place_holder_size=10):
        paths = []
        stack = [[x.contig_id] for x in self.contigs]
        best_costs = [float("Inf")] * place_holder_size  # five best costs
        while stack:
            best_costs.sort()
            path = stack.pop()
            cost = self.compute_assembly_cost([self.get_contig(x) for x in path])
            # Trimming conditions
            if cost >= best_costs[-1] and not cost == float("Inf"):
                # cost can only get worst, so trim this search
                continue

            # if path > query_length

            n = path[-1]
            if n in self.assembly_graph:
                for node in self.assembly_graph[n]:
                    new_path = path[:] + [node]
                    stack.append(new_path)
            else:
                if cost < best_costs[-1]:
                    best_costs[-1] = cost
                paths.append(path)
        return paths
