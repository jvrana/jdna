'''
Project: jdna
File: das_assembly
Author: Justin
Date: 2/13/17

Description: 

'''

from das_contig import *


class AssemblyGraph(ContigContainer):
    def __init__(self, primers=None, contigs=None):
        super(AssemblyGraph, self).__init__(meta=contigs.meta.__dict__, contigs=contigs.contigs)
        self.graph = {}
        self.primers = primers
        if primers is None:
            self.primers = ContigContainer()
        self.make_dictionary()
        self.assemblies = []

    def get_all_assemblies(self, sort=True, place_holder_size=5, save_history=False):
        self.make_assembly_graph(sort=sort)
        self.assemblies = None
        if save_history:
            self.assemblies, self.assembly_history = self.dfs_iter(place_holder_size=place_holder_size, data_plot=True)
        else:
            self.assemblies = self.dfs_iter(place_holder_size=place_holder_size, data_plot=False)
        self.assemblies = sorted(self.assemblies, key=lambda x: x.total_cost())
        return self.assemblies

    def make_assembly_graph(self, sort=True):
        self.make_dictionary()
        graph = {}
        pairs = list(itertools.permutations(self.contigs, 2))
        print '{} pairs to search'.format(len(pairs))
        num_nodes = 0
        for l, r in pairs:

            if l == r:
                continue
            if Assembly.assembly_condition(l, r):
                if l.contig_id not in graph:
                    graph[l.contig_id] = []
                graph[l.contig_id].append(r.contig_id)
                num_nodes += 1
        if sort:
            for k in graph.keys():
                a_array = graph[k][:]

                # graph[k] = sorted(a_array, key=lambda x: self.compute_assembly_cost(
                #     [self.get_contig(k), self.get_contig(x)])[0])

                graph[k] = sorted(a_array, key=lambda x: self.get_contig(x).query.end - self.get_contig(x).query.start,
                                  reverse=True)
        self.graph = graph
        print 'Graph Size: {}'.format(num_nodes)
        # self.assembly_graph['root'] = [x.contig_id for x in self.contigs]
        return self.graph

    def dfs_iter(self, place_holder_size=5, data_plot=False):
        assemblies = []
        best_costs_array = []
        steps = []
        step = 0
        sorted_contigs = sorted(self.contigs, key=lambda x: x.query.end - x.query.start, reverse=True)
        stack = []
        for c in sorted_contigs:
            a = Assembly([c], self, self.primers)
            stack.append(a)
        best_costs = [float("Inf")] * place_holder_size  # five best costs
        while stack:
            step += 1
            steps.append(step)
            best_costs.sort()
            assembly = stack.pop()

            gap = assembly.get_gaps()
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
                if data_plot:
                    best_costs_array.append(best_costs[:])
                print best_costs
            assembly.best_cost = assembly.total_cost()
            assemblies.append(assembly)
            # end of path
        assemblies = sorted(assemblies,
                            key=lambda x: x.total_cost())
        if data_plot:
            return assemblies, zip(steps, best_costs_array)
        return assemblies

    def dump(self, out):
        j = deepcopy(self.__dict__)
        del j['contig_dictionary']
        j['meta'] = j['meta'].__dict__
        j['contigs'] = [x.json() for x in j['contigs']]
        del j['meta']['query_seq']
        del j['meta']['contig_seqs']
        del j['primers']
        with open(out, 'w') as output:
            json.dump(j, output)

# TODO: share parameters from a single file
class Assembly(ContigContainer):
    MIN_PCR_SIZE = 250.
    MAX_PCR_SIZE = 5000.
    MAX_HOMOLOGY = 100.
    MIN_ASSEMBLY_HOMOLOGY = 20.0 # TODO: implement this parameter
    PRIMER_COST = 15.
    PCR_COST = 14.
    FIVEPRIME_EXT_REACH = 20.  # 'reachability' for a extended primer
    SYNTHESIS_THRESHOLD = 125.  # min bp for synthesis
    SYNTHESIS_COST = 0.11  # per bp
    SYNTHESIS_MIN_COST = 89.  # per synthesis
    assembly_id = 0

    def __init__(self, assembly_path, contig_container, primer_container):
        super(Assembly, self).__init__(meta=contig_container.meta.__dict__, contigs=None)
        self.id = Assembly.assembly_id
        Assembly.assembly_id += 1
        self.contig_container = contig_container
        self.primer_container = primer_container
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
        return self.assembly_span() < self.meta.query_length + Assembly.MAX_HOMOLOGY

    @staticmethod
    def gap(left, right):
        """
        Calculates the unreachable gap
        :param left:
        :param right:
        :return:
        """
        d = right.query.start - left.query.end
        r = 0
        if right.end_label == Contig.NEW_PRIMER:
            r += Assembly.FIVEPRIME_EXT_REACH
        if left.start_label == Contig.NEW_PRIMER:
            r += Assembly.FIVEPRIME_EXT_REACH
        gap = d - r
        return gap

    def circular_gap(self):
        """
        Calculates the gap for a circular contig path.
        Note: non-increasing gap function; cannot be used for graph trimming
        :return: calculated gap between the last and first contig for circularization
        """

        return self.meta.query_length + Assembly.gap(self.contigs[-1], self.contigs[0])

    def print_contig_path(self):
        c_arr = []
        for c in self.contigs:
            if isinstance(c, int):
                c = self.get_contig(c)
            c_arr.append(c)

        print [(c.query.start, c.query.end) for c in c_arr]

    def get_assembly_pairs(self):
        if len(self.contigs) == 1:
            return []
        pairs = zip(self.contigs[:-1], self.contigs[1:])
        return pairs

    def get_gaps(self):
        """
        Calculates all non-circular gaps for a contig_path
        :return:
        """
        gaps = []
        for l, r in self.get_assembly_pairs():
            if Assembly.assembly_condition(l, r):
                gaps.append(Assembly.gap(l, r))
            else:
                gaps.append(float("Inf"))
                return gaps
        return gaps

    # TODO: make this more complicated and reflective of the actual gap cost
    # g = [
    #         (0, [-Assembly.MAX_HOMOLOGY, -Assembly.MIN_ASSEMBLY_HOMOLOGY]),
    #     (PRIMER_COST, [-Assembly.MIN_ASSEMBLY_HOMOLOGY+1, -10]),
    #     (2X Priemrs),
    #     (LONG PRIMERS),
    #     etc.
    #
    def _gap_cost(self, gap):
        if gap <= -Assembly.MIN_ASSEMBLY_HOMOLOGY:
            return 0
        else:
            return (gap + Assembly.MIN_ASSEMBLY_HOMOLOGY) * Assembly.SYNTHESIS_COST

    def get_gap_cost(self):
        """
        A monotonically increasing gap cost function.
        :return:
        """
        gaps = self.get_gaps()
        return sum([self._gap_cost(g) for g in gaps])

    @staticmethod
    def get_pcr_cost(contig):
        """
        Cost of a pcr
        :param contig:
        :return:
        """

        if not Assembly.MIN_PCR_SIZE < contig.query.end - contig.query.start < Assembly.MAX_PCR_SIZE:
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
        if cost < 0:
            cost = 0
        return cost

    def assembly_span(self):
        return self.last().query.end - self.first().query.start

    def unassembled_span(self):
        return self.meta.query_length - self.assembly_span()

    def unassembled_span_cost(self):
        return self._gap_cost(self.circular_gap())

    def get_num_synthesis_fragments(self):
        gaps = self.get_gaps()
        gaps.append(self.circular_gap())
        return len(filter(lambda x: x > Assembly.SYNTHESIS_THRESHOLD, gaps))

    def new_synthesis_cost(self):
        return self.get_num_synthesis_fragments() * Assembly.SYNTHESIS_MIN_COST + \
               self.get_gap_cost() + self.unassembled_span_cost()

    def get_fragment_cost(self):
        return sum([Assembly.get_pcr_cost(c) for c in self.contigs])

    @property
    def assembly_probability(self):
        """
        The probability of a successful GIBSON/SLIC assembly given
        a number of fragments.
        """

        K = 6.0
        n = 3.0

        p = 1/(1+(len(self.contigs)/K)**n)
        if p <= 0:
            p = 0.0001
        return p

    def total_cost(self):
        """
        Calculates total cost in dollars for a given assembly.
        :return:
        """
        self.cost = (self.new_synthesis_cost() + self.get_fragment_cost()) / self.assembly_probability
        return self.cost

    @staticmethod
    def assembly_condition(left, right):
        """
        Determines whether two contigs can be assembled either directly or
        with newly syntheisized intermediates. Used loosely to assemble the contig
        graph.
        :param left:
        :param right:
        :return:
        """
        r_5prime_threshold = 0
        r_3prime_threshold = Assembly.MAX_HOMOLOGY
        l_pos = left.query.end
        r_pos = right.query.start
        # return r_pos == l_pos + 1 # if its consecutive
        return r_pos > l_pos - r_3prime_threshold and \
               right.query.end > left.query.end

    def get_all_templates(self):
        filenames = []
        for c in self.contigs:
            filenames.append(c.filename)
        return list(set(filenames))

    def dump(self, out):
        j = deepcopy(self.__dict__)
        del j['contig_dictionary']
        j['meta'] = j['meta'].__dict__
        j['contigs'] = [x.json() for x in j['contigs']]
        del j['meta']['query_seq']
        del j['meta']['contig_seqs']
        del j['contig_container']
        del j['primer_container']
        with open(out, 'w') as output:
            json.dump(j, output)

    def summary(self):
        contigs = ''
        for c in self.contigs:
            contigs += '''
        \tContig {id}
        \t\tCost: {cost}
        \t\tQuery Span: {start}-{end}
        \t\tRel Loc: {rs}-{re}
        \t\tSub Loc: {s_start}-{s_end}
        \t\tSubject Acc: {subject}
        \t\tPrimers: {fwd}, {rev}
        '''.format(id=c.contig_id,
                   cost=self.get_pcr_cost(c),
                   start=c.query.start,
                   end=c.query.end,
                   rs=c.query.start - self.contigs[0].query.start,
                   re=c.query.end - self.contigs[0].query.start,
                   subject=c.subject.name,
                   fwd=c.start_label,
                   rev=c.end_label,
                   s_start=c.subject.start,
                   s_end=c.subject.end
                   )

        summary_str = '''
** Assembly Summary **
ID: {id}
Query: {query}
Query Length: {querylength}
Contig Breakdown
\tNum Frags: {num}
{contigs}
Breakdown {totalcost}
\tGap Cost: {gaps}, ${gapcost}
\tAssembled Span: {span}
\tUnassembled Span: {uspan}, ${uspancost}
\tFragment Costs: ${fragcost}, {fragbreakdown}
\tNew Synthesis Costs: ${newsynth}
\tProbability: {probability}

        '''.format(
            id=self.id,
            num=len(self.contigs),
            totalcost=self.total_cost(),
            gaps=self.get_gaps(),
            gapcost=self.get_gap_cost(),
            contigs=contigs,
            span=self.assembly_span(),
            uspan=self.circular_gap(),
            uspancost=self.unassembled_span_cost(),
            fragcost=self.get_fragment_cost(),
            fragbreakdown=[self.get_pcr_cost(c) for c in self.contigs],
            query=self.meta.query,
            querylength=self.meta.query_length,
            newsynth=self.new_synthesis_cost(),
            probability=self.assembly_probability
        )
        return summary_str