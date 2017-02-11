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


class QueryRegion(object):
    def __init__(self, **kwargs):
        self.query_acc = kwargs['query_acc']
        self.q_start = kwargs['q_start']
        self.q_end = kwargs['q_end']
        self.query_length = kwargs['query_length']

    def get_span(self):
        return self.q_end - self.q_start


    def json(self):
        return self.__dict__

class Contig(QueryRegion):
    contig_id = 0
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
        '''

        :param startpoints: list of (pos, label) tuples
        :param endpoints: list of (pos, label) tuples
        :param include_contig: whether to include the contig endpoints
        :param contig_type: contig_label, used to determine the type of contig this is for the cost analysis
        :param circular: used to determine the type of contig this is for the cost analysis
        :return:
        '''
        positions = list(itertools.product(startpoints, endpoints))
        contigs = []
        for start, end in positions:
            if isinstance(start, int):
                x = 1
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

    # TODO: add primer label names to start_label and end_label
    def pcr_products_of_contig(self, primers):
        '''
        Produces all possible fragments from primer positions
        :param contig:
        :param primers:
        :param ignore_direction:
        :return:
        '''
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
            minimum_primer_anneal = 15  # TODO: move MIN_PRIMER to design parameters
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
        j['contigs'] = [x.json() for x in j['contigs']]
        del j['meta']['query_seq']
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


class AssemblyGraph(ContigContainer):
    def __init__(self, primers=None, contigs=None):
        super(AssemblyGraph, self).__init__(meta=contigs.meta.__dict__, contigs=contigs.contigs)
        self.graph = {}
        self.primers = primers
        if primers is None:
            self.primers = ContigContainer()
        self.make_dictionary()

    def get_all_assemblies(self, sort=True, place_holder_size=5):
        self.make_assembly_graph(sort=sort)
        self.assemblies = self.dfs_iter(place_holder_size=place_holder_size)
        self.assemblies = sorted(self.assemblies, key=lambda x: x.total_cost())
        return self.assemblies

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

                graph[k] = sorted(a_array, key=lambda x: self.get_contig(x).q_end - self.get_contig(x).q_start,
                                  reverse=True)
        self.graph = graph
        # self.assembly_graph['root'] = [x.contig_id for x in self.contigs]
        return self.graph

    def dfs_iter(self, place_holder_size=5):
        assemblies = []
        sorted_contigs = sorted(self.contigs, key=lambda x: x.q_end - x.q_start, reverse=True)
        stack = []
        for c in sorted_contigs:
            a = Assembly([c], self, self.primers)
            stack.append(a)
        best_costs = [float("Inf")] * place_holder_size  # five best costs
        while stack:
            best_costs.sort()
            assembly = stack.pop()

            gap_cost = assembly.get_gap_cost()
            frag_cost = assembly.get_fragment_cost()

            cost = assembly.total_cost()

            # Trimming conditions
            # Gap and frag cost are monotonically increasing and so are used for trimming.
            if gap_cost + frag_cost > best_costs[-1] and not cost == float("Inf"):
                # cost can only get worst, so trim this search
                continue
            if cost == float("Inf"):
                continue

            can_extend = False
            has_children = False

            new_paths = []
            n = assembly.last().contig_id
            if n in self.graph:
                children = self.graph[n]
                if children:
                    has_children = True
                    for node in children:
                        assembly_copy = copy(assembly)
                        assembly_copy.contigs = assembly.contigs[:] + [self.get_contig(node)]
                        new_paths.append(assembly_copy)

            for new_path in new_paths:
                if new_path.can_extend():
                    can_extend = True
                    stack.append(new_path)
            if cost < best_costs[-1]:
                best_costs[-1] = cost
            assemblies.append(assembly)
            # end of path
        assemblies = sorted(assemblies,
                            key=lambda x: x.total_cost())
        return assemblies


class Assembly(ContigContainer):
    MIN_PCR_SIZE = 250.
    MAX_PCR_SIZE = 10000.
    PRIMER_COST = 15.
    PCR_COST = 14.
    FIVEPRIME_EXT_REACH = 20.  # 'reachability' for a extended primer
    SYNTHESIS_THRESHOLD = 0.  # min bp for synthesis
    SYNTHESIS_COST = 0.11  # per bp
    SYNTHESIS_MIN_COST = 89.  # per synthesis

    def __init__(self, assembly_path, contig_container, primer_container, meta=None):
        super(Assembly, self).__init__(meta=contig_container.meta.__dict__, contigs=None)
        self.contig_container = contig_container
        self.primer_continer = primer_container
        self._convert_assembly_path(assembly_path)

    def _convert_assembly_path(self, assembly_path):
        self.contig_container.make_dictionary()
        for i, ap in enumerate(assembly_path):
            if isinstance(ap, Contig):
                pass
            elif isinstance(ap, int):
                assembly_path[i] = self.contig_container.get_contig(ap)
        self.contigs = assembly_path

    def first(self):
        return self.contigs[0]

    def last(self):
        return self.contigs[-1]

    def can_extend(self):
        return self.assembly_span() < self.meta.query_length

    def fill_contig_gaps(self):
        pairs = self.get_assembly_pairs()
        new_contigs = []
        for l, r in pairs:
            new_contigs.append(l)
            q = QueryRegion(query_acc=l.query_acc, query_length=l.query_length, q_start=l.q_end, q_end=r.q_start)
            if q.get_span() > 0:
                new_contigs.append(q)
        new_contigs.append(pairs[-1][1])  # append last contig
        self.contigs = new_contigs

    @staticmethod
    def gap(left, right):
        '''
        Calculates the unreachable gap
        :param left:
        :param right:
        :return:
        '''
        d = right.q_start - left.q_end
        r = 0
        if right.end_label == Contig.NEW_PRIMER or right.end_label is Contig.DIRECT_END:
            r += Assembly.FIVEPRIME_EXT_REACH
        if left.start_label == Contig.NEW_PRIMER or left.start_label is Contig.DIRECT_END:
            r += Assembly.FIVEPRIME_EXT_REACH
        gap = d - r
        return gap

    @staticmethod
    def circular_gap(left, right, query_length):
        '''
        Calculates the gap for a circular contig path.
        Note: non-increasing gap function; cannot be used for graph trimming
        :param left:
        :param right:
        :param query_length:
        :return:
        '''
        return query_length + Assembly.gap(right, left)

    def print_contig_path(self):
        c_arr = []
        for c in self.contigs:
            if isinstance(c, int):
                c = self.get_contig(c)
            c_arr.append(c)

        print [(c.q_start, c.q_end) for c in c_arr]

    def get_assembly_pairs(self):
        return zip(self.contigs[:-1], self.contigs[1:])

    def get_gaps(self):
        '''
        Calculates all non-circular gaps for a contig_path
        :param contig_path:
        :param circular:
        :return:
        '''
        gaps = []
        for l, r in self.get_assembly_pairs():
            if Assembly.assembly_condition(l, r):
                gaps.append(Assembly.gap(l, r))
            else:
                gaps.append(float("Inf"))
                return gaps
        return gaps

    def get_gap_cost(self):
        '''
        A monotonically increasing gap cost function.
        :param gaps:
        :return:
        '''
        gaps = self.get_gaps()
        non_zero_gaps = filter(lambda x: x > 0, gaps)
        return sum(non_zero_gaps) * Assembly.SYNTHESIS_COST

    @staticmethod
    def pcr_cost(contig):
        '''
        Cost of a pcr
        :param contig:
        :return:
        '''

        if not Assembly.MIN_PCR_SIZE < contig.q_end - contig.q_start < Assembly.MAX_PCR_SIZE:
            return float("Inf")
        if contig.is_direct():
            # direct_synthesis
            return 0
        cost = Assembly.PCR_COST
        if contig.circular:
            cost += 2 * Assembly.PRIMER_COST
        else:
            if contig.start_label == Contig.NEW_PRIMER or contig.start_label is Contig.DIRECT_END:
                cost += Assembly.PRIMER_COST
            if contig.end_label == Contig.NEW_PRIMER or contig.end_label is Contig.DIRECT_END:
                cost += Assembly.PRIMER_COST
        # print contig.circular, contig.parent, contig.start_label, contig.end_label, cost
        return cost

    def assembly_span(self):
        return self.last().q_end - self.first().q_start

    def unassembled_span(self):
        return self.meta.query_length - self.assembly_span()

    def unassembled_span_cost(self):
        return self.unassembled_span() * Assembly.SYNTHESIS_COST

    def get_num_synthesis_fragments(self):
        gaps = self.get_gaps()
        gaps.append(self.unassembled_span())
        return len(filter(lambda x: x > Assembly.SYNTHESIS_THRESHOLD, gaps))

    def new_synthesis_cost(self):
        return self.get_num_synthesis_fragments() * Assembly.SYNTHESIS_MIN_COST + \
               self.get_gap_cost() + self.unassembled_span_cost()

    def get_fragment_cost(self):
        return sum([Assembly.pcr_cost(c) for c in self.contigs])

    def assembly_probability(self):
        '''
        The probability of a successful GIBSON/SLIC assembly given
        a number of fragments.
        :param num_frags:
        :return:
        '''
        p = 1.0 - (1.0 / 8.0) * len(self.contigs)
        if p <= 0:
            p = 0.0001
        return p

    def total_cost(self):
        '''
        Calculates total cost in dollars for a given assembly.
        :param contig_path:
        :param query_length:
        :param circular:
        :return:
        '''
        return (self.new_synthesis_cost() + self.get_fragment_cost()) / self.assembly_probability()

    @staticmethod
    def assembly_condition(left, right):
        '''
        Determines whether two contigs can be assembled either directly or
        with newly syntheisized intermediates. Used loosely to assemble the contig
        graph.
        :param left:
        :param right:
        :return:
        '''
        r_5prime_threshold = 0
        r_3prime_threshold = 60
        l_pos = left.q_end
        r_pos = right.q_start
        # return r_pos == l_pos + 1 # if its consecutive
        return r_pos > l_pos - r_3prime_threshold and \
               right.q_end > left.q_end

    def get_all_templates(self):
        filenames = []
        for c in self.contigs:
            filenames.append(c.filename)
        return list(set(filenames))


class J5Assembly(Assembly):
    PARTLABELS = ['Part Source (Sequence Display ID)', 'Part Name', 'Reverse Compliment?', 'Start (bp)', 'End (bp)',
                  'Five Prime Internal Preferred Overhangs?', 'Three Prime Internal Preferred Overhangs?']
    SEQLISTLABELS = ['Sequence File Name', 'Format']
    TARGETLABELS = ["(>Bin) or Part Name", "Direction", "Forced Assembly Strategy?",
                    "Forced Relative Overhang Position?", "Direct Synthesis Firewall?", "Extra 5' CPEC overlap bps",
                    "Extra 3' CPEC overlap bps"]
    MASTEROLIGOSLABELS = ['Oligo Name', 'Length', "Tm", "Tm (3' only)", "Sequence"]
    MASTERPLASMIDLABELS = ["Plasmid Name", "Alias", "Contents", "Length", "Sequence"]
    MASTERDIRECTLABELS = ["Direct Synthesis Name", "Alias", "Contents", "length", "Sequence"]
    ZIP = None
    EUGENE = None
    params = None

    def __init__(self, assembly):
        super(J5Assembly, self).__init__(assembly.contigs, assembly.contig_container, assembly.primer_continer)

    def to_csv(self, labels, rows):
        a = [labels] + rows
        csv = '\n'.join([str(x) for x in a])
        return csv

    def encode64(self, filestring):
        return base64.encodestring(filestring)

    def target(self):
        rows = []
        for c in self.contigs:
            row = [
                c.contig_id,
                'forward',
                '',
                '',
                '',
                '',
                ''
            ]
            rows.append(row)
        self.target = self.encode64(self.to_csv(J5Assembly.TARGETLABELS, rows))

    def parts(self):
        rows = []
        for c in self.contigs:
            row = [
                c.contig_id,
                c.seqrecord.id,
                str(False).upper(),
                c.s_start,
                c.s_end,
                '',
                ''
            ]
            rows.append(row)
        csv = self.to_csv(J5Assembly.PARTLABELS, rows)
        self.parts = self.encode64(csv)

    def direct(self):
        self.direct = self.encode64(self.to_csv(J5Assembly.MASTERDIRECTLABELS, []))

    def plasmids(self):
        self.plasmids = self.encode64(self.to_csv(J5Assembly.MASTERPLASMIDLABELS, []))

    def primers(self):
        self.primers = self.encode64(self.to_csv(J5Assembly.MASTEROLIGOSLABELS, []))

    def sequences(self):
        sequences = []
        for c in self.contigs:
            if c.is_direct():
                pass
            sequences.append((c.seqrecord.id, c.filename))
        rows = list(set(sequences))

        # Save sequence_list
        csv = self.to_csv(J5Assembly.SEQLISTLABELS, rows)
        self.seq_list = self.encode64(csv)

        # TODO: move zip stuff to seqio
        # Save zipped sequence file
        temp = tempfile.TemporaryFile() #open('/Users/Justin/Desktop/temporary.txt', 'w')#tempfile.TemporaryFile()
        with zipfile.ZipFile(temp, 'w') as zf:

            for s in sequences:
                filename = s[1]
                with open(filename, 'rU') as handle:
                    zf.writestr(os.path.join('ZippedTemplates', os.path.basename(filename)), handle.read())
        temp.seek(0)
        self.encoded_zip = self.encode64(temp.read())
        temp.close()
