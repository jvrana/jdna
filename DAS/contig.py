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
        self.__dict__.update(kwargs)


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
        self.parent = 'query'
        self.end_label = 'query'
        self.start_label = 'query'
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

    def divide_contig(self, startpoints, endpoints, include_contig=True, contig_type=None, circular=None, start_label='', end_label=''):
        positions = list(itertools.product(startpoints, endpoints))
        contigs = []
        for s, e in positions:
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
                    new_contig.start_label = 'contig'
                else:
                    new_contig.start_label = start_label
                if new_contig.q_end == self.q_end:
                    new_contig.end_label = 'contig'
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

    def is_perfect(self):
        return self.alignment_length == self.identical and \
                self.gaps == 0 and self.gap_opens == 0

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

    def get_contig(self, id):
        return self.contig_dictionary[id]

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
        j['contigs'] = [x.__dict__ for x in j['contigs']]
        with open(out, 'w') as output:
            json.dump(j, output)

    def sort_contigs(self):
        self.contigs = sorted(self.contigs, key=lambda x: x.q_start)

    def filter_perfect(self):
        print len(self.contigs),
        filtered_contigs = []
        for contig in self.contigs:
            if contig.is_perfect():
                filtered_contigs.append(contig)
            else:
                print contig.__dict__
        print len(self.contigs)

    # TODO: handle ambiquoous NNNN dna in blast search by eliminating gap_opens, gaps if they are N's


class Assembly(ContigContainer):

    def __init__(self, meta=None, contigs=[]):
        super(Assembly, self).__init__(meta=meta, contigs=contigs)
        self.assembly_graph = {}
        self.paths = []

    def get_all_assemblies(self, sort=True, place_holder_size=5):
        self.make_assembly_graph(sort=sort)
        self.assemblies = self.dfs_iter(place_holder_size=place_holder_size)
        return self.assemblies

    def print_contig_path(self, contigs):
        c_arr = []
        for c in contigs:
            if isinstance(c, int):
                c = self.get_contig(c)
            c_arr.append(c)

        print [(c.q_start, c.q_end) for c in c_arr]



    @staticmethod
    def assembly_condition(left, right):
        r_5prime_threshold = 0
        r_3prime_threshold = 60
        l_pos = left.q_end
        r_pos = right.q_start
        # return r_pos == l_pos + 1 # if its consecutive
        return r_pos > l_pos - r_3prime_threshold and \
               right.q_end > left.q_end

    def make_assembly_graph(self, sort=True):
        self.make_dictionary()
        graph = {}
        pairs = list(itertools.permutations(self.contigs, 2))
        print '{} pairs to search'.format(len(pairs))
        for l, r in pairs:

            if l == r:
                continue
            if Assembly.assembly_condition(l, r):
                if l.contig_id not in graph:
                    graph[l.contig_id] = []
                graph[l.contig_id].append(r.contig_id)

        if sort:
            for k in graph.keys():
                a_array = graph[k][:]

                # graph[k] = sorted(a_array, key=lambda x: self.compute_assembly_cost(
                #     [self.get_contig(k), self.get_contig(x)])[0])

                graph[k] = sorted(a_array, key=lambda x: self.get_contig(x).q_end - self.get_contig(x).q_start, reverse=True)
        self.assembly_graph = graph
        # self.assembly_graph['root'] = [x.contig_id for x in self.contigs]
        return self.assembly_graph

    @staticmethod
    def fragment_cost(contig):
        MIN_PCR_SIZE = 250
        MAX_PCR_SIZE = 10000
        PRIMER_COST = 15
        PCR_COST = 14

        if not MIN_PCR_SIZE < contig.q_end - contig.q_start < MAX_PCR_SIZE:
            return float("Inf")
        if not contig.circular and contig.parent == 'query':
            # direct_synthesis
            return 0
        cost = 15
        if contig.circular:
            cost += 2*PRIMER_COST
        else:
            if not contig.start_label == 'primer':
                cost += PRIMER_COST
            if not contig.end_label == 'primer':
                cost += PRIMER_COST
        # print contig.circular, contig.parent, contig.start_label, contig.end_label, cost
        return cost

    def compute_assembly_cost(self, contig_path):
        if len(contig_path) == 0:
            return float("Inf")

        pairs = zip(contig_path[:-1], contig_path[1:])


        query_span = contig_path[-1].q_end - contig_path[0].q_start
        span_cost = self.meta.query_length - query_span
        if span_cost < 0:
            span_cost = float("Inf")
        gaps = []
        num_synthesized_fragments = 0
        for l, r in pairs:
            if Assembly.assembly_condition(l, r): # valid assembly
                d = r.q_start - l.q_end
                reachability = 0
                if r.end_label == 'new primer':
                    reachability += 20
                if l.end_label == 'new primer':
                    reachability += 20
                gap_distance = d - reachability
                if gap_distance < 0:
                    gap_distance = 0
                if gap_distance > 60:
                    num_synthesized_fragments += 1
                gaps.append(gap_distance)
            else: # non-valid assembly
                gaps += [float("Inf")]
                break
                # non valid assembly

        non_zero_gaps = filter(lambda x: x != 0, gaps)
        # non_zero_gaps = filter(lambda x: x > 60, gaps)
        gap_cost = sum(non_zero_gaps)
        if span_cost > 0:
            r, l = contig_path[-1], contig_path[0]
            d = r.q_start - l.q_end
            reachability = 0
            if r.end_label == 'new primer':
                reachability += 20
            if l.end_label == 'new primer':
                reachability += 20
            gap_distance = d - reachability
            if gap_distance < 0:
                gap_distance = 0
            if gap_distance > 60:
                num_synthesized_fragments += 1
        fragment_cost = sum([Assembly.fragment_cost(c) for c in contig_path])
        return gap_cost, span_cost, num_synthesized_fragments, fragment_cost

    def total_cost(self, contig_path):
        gap_cost, span_cost, new_frags, fragment_cost = self.compute_assembly_cost(contig_path)
        probability = 1.0 - (1.0/8.0) * (new_frags + len(contig_path))
        if probability <= 0:
            probability = 0.001
        return (1/(probability) *(gap_cost + span_cost) * 0.11 + fragment_cost + 89*new_frags)

    def dfs_iter(self, place_holder_size=5):
        paths = []
        sorted_contigs = sorted(self.contigs, key=lambda x: x.q_end - x.q_start, reverse=True)
        stack = [[x.contig_id] for x in sorted_contigs]
        best_costs = [float("Inf")] * place_holder_size  # five best costs
        while stack:
            best_costs.sort()
            path = stack.pop()
            contig_path = [self.get_contig(x) for x in path]
            gap_cost, span_cost, num_frags, frag_cost = self.compute_assembly_cost(contig_path)
            cost = self.total_cost(contig_path)
            # Trimming conditions
            if gap_cost*0.11 + frag_cost > best_costs[-1] and not cost == float("Inf"):
                # cost can only get worst, so trim this search
                continue
            if cost == float("Inf"):
                continue

            # if path > query_length

            n = path[-1]
            can_extend = False
            has_children = False

            new_paths = []
            if n in self.assembly_graph:
                children = self.assembly_graph[n]
                if children:
                    has_children = True
                    for node in children:
                        new_paths.append(path[:] + [node])

            for new_path in new_paths:
                query_span = self.get_contig(new_path[-1]).q_end - self.get_contig(new_path[0]).q_start
                if query_span < self.meta.query_length:
                    can_extend = True
                    stack.append(new_path)

            if not (has_children and can_extend):
                if cost < best_costs[-1]:
                    best_costs[-1] = cost
                    print cost, best_costs
                    paths.append(path)
            # end of path
        paths = sorted(paths, key=lambda x: self.total_cost([self.get_contig(c) for c in x]))
        return paths

    def contigs_query_span(self, contig_list):
        contigs = []
        for c in contigs:
            if isinstance(c, int):
                c = self.get_contig(c)
            contigs.append(c)

        sort_by_start = sorted(contigs, key=lambda x: x.q_start)
        sort_by_end = sorted(contigs, key=lambda x: x.q_end)
        return sort_by_end[-1].q_end - sort_by_start[0].q_start

    @staticmethod
    def get_primers_within_bounds(contig, primers, minimum_primer_anneal=15):
        def primer_within_bounds(primer):
            minimum_primer_anneal = 15  # TODO: move MIN_PRIMER to design parameters
            if primer.subject_strand == 'plus':
                return contig.q_start + minimum_primer_anneal < primer.q_end < contig.q_end
            elif primer.subject_strand == 'minus':
                return contig.q_start < primer.q_start < contig.q_end - minimum_primer_anneal
            return False

        return filter(lambda x: primer_within_bounds(x), primers)

    @staticmethod
    def pcr_products_of_contig(contig, primers):
        '''
        Produces all possible fragments from primer positions
        :param contig:
        :param primers:
        :param ignore_direction:
        :return:
        '''

        new_alignments = []
        primers = Assembly.get_primers_within_bounds(contig, primers)
        rev_primer_pos = []
        fwd_primer_pos = []
        new_fwd_primer_pos = []
        new_rev_primer_pos = []
        for primer in primers:
            direction = primer.subject_strand
            if direction == 'plus':
                fwd_primer_pos.append(primer.q_start)
                new_rev_primer_pos.append(primer.q_end)
            if direction == 'minus':
                rev_primer_pos.append(primer.q_end)
                new_fwd_primer_pos.append(primer.q_start)
        contigs = []
        contigs += contig.divide_contig(fwd_primer_pos, rev_primer_pos, include_contig=False,
                                    contig_type='product', circular=False, start_label='primer', end_label='primer')
        contigs += contig.divide_contig(fwd_primer_pos, new_rev_primer_pos + [contig.q_end], include_contig=False,
                                        contig_type='product', circular=False, start_label='primer', end_label='new_primer')
        contigs += contig.divide_contig(new_fwd_primer_pos + [contig.q_start], rev_primer_pos, include_contig=False,
                                        contig_type='product', circular=False, start_label='new_primer',
                                        end_label='primer')
        contigs += contig.divide_contig(new_fwd_primer_pos + [contig.q_start], new_rev_primer_pos + [contig.q_end], include_contig=False,
                                        contig_type='product', circular=False, start_label='new_primer',
                                        end_label='new_primer')
        return contigs


    # TODO: generate additional pcr products that could be homologous to existing primer
    # TODO: compute alignment_graph for each 'contig'

    def expand_contigs(self, primers):
        all_alignments = []
        for contig in self.contigs:
            new_alignments = Assembly.pcr_products_of_contig(contig, primers)
            all_alignments += new_alignments
        self.contigs += all_alignments