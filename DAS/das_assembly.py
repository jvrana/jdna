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
            assemblies.append(assembly)
            # end of path
        assemblies = sorted(assemblies,
                            key=lambda x: x.total_cost())
        return assemblies

# TODO: share parameters from a single file
class Assembly(ContigContainer):
    MIN_PCR_SIZE = 250.
    MAX_PCR_SIZE = 10000.
    PRIMER_COST = 15.
    PCR_COST = 14.
    FIVEPRIME_EXT_REACH = 20.  # 'reachability' for a extended primer
    SYNTHESIS_THRESHOLD = 0.  # min bp for synthesis
    SYNTHESIS_COST = 0.11  # per bp
    SYNTHESIS_MIN_COST = 89.  # per synthesis
    assembly_id = 0

    def __init__(self, assembly_path, contig_container, primer_container):
        super(Assembly, self).__init__(meta=contig_container.meta.__dict__, contigs=None)
        self.id = Assembly.assembly_id
        Assembly.assembly_id += 1
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
                print q
                new_contigs.append(q)
        new_contigs.append(pairs[-1][1])  # append last contig
        self.contigs = new_contigs

    @staticmethod
    def gap(left, right):
        """
        Calculates the unreachable gap
        :param left:
        :param right:
        :return:
        """
        d = right.q_start - left.q_end
        r = 0
        if right.end_label == Contig.NEW_PRIMER or right.end_label is Contig.DIRECT_END:
            r += Assembly.FIVEPRIME_EXT_REACH
        if left.start_label == Contig.NEW_PRIMER or left.start_label is Contig.DIRECT_END:
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

        print [(c.q_start, c.q_end) for c in c_arr]

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

    def get_gap_cost(self):
        """
        A monotonically increasing gap cost function.
        :return:
        """
        gaps = self.get_gaps()
        non_zero_gaps = filter(lambda x: x > 0, gaps)
        return sum(non_zero_gaps) * Assembly.SYNTHESIS_COST

    @staticmethod
    def pcr_cost(contig):
        """
        Cost of a pcr
        :param contig:
        :return:
        """

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
        gaps.append(self.circular_gap())
        return len(filter(lambda x: x > Assembly.SYNTHESIS_THRESHOLD, gaps))

    def new_synthesis_cost(self):
        return self.get_num_synthesis_fragments() * Assembly.SYNTHESIS_MIN_COST + \
               self.get_gap_cost() + self.unassembled_span_cost()

    def get_fragment_cost(self):
        return sum([Assembly.pcr_cost(c) for c in self.contigs])

    @property
    def assembly_probability(self):
        """
        The probability of a successful GIBSON/SLIC assembly given
        a number of fragments.
        """
        p = 1.0 - (1.0 / 8.0) * len(self.contigs)
        if p <= 0:
            p = 0.0001
        return p

    def total_cost(self):
        """
        Calculates total cost in dollars for a given assembly.
        :return:
        """
        return (self.new_synthesis_cost() + self.get_fragment_cost()) / self.assembly_probability

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
                   cost=self.pcr_cost(c),
                   start=c.q_start,
                   end=c.q_end,
                   rs=c.q_start - self.contigs[0].q_start,
                   re=c.q_end - self.contigs[0].q_start,
                   subject=c.subject_acc,
                   fwd=c.start_label,
                   rev=c.end_label,
                   s_start=c.s_start,
                   s_end=c.s_end
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
\t

        '''.format(
            id=self.id,
            num=len(self.contigs),
            totalcost=self.total_cost(),
            gaps=self.get_gaps(),
            gapcost=self.get_gap_cost(),
            contigs=contigs,
            span=self.assembly_span(),
            uspan=self.unassembled_span(),
            uspancost=self.unassembled_span_cost(),
            fragcost=self.get_fragment_cost(),
            fragbreakdown=[self.pcr_cost(c) for c in self.contigs],
            query=self.meta.query,
            querylength=self.meta.query_length,
            newsynth=self.new_synthesis_cost()
        )
        return summary_str


class J5Assembly(Assembly):
    PARTLABELS = ['Part Name', 'Part Source (Sequence Display ID)', 'Reverse Compliment?', 'Start (bp)', 'End (bp)',
                  'Five Prime Internal Preferred Overhangs?', 'Three Prime Internal Preferred Overhangs?']
    SEQLISTLABELS = ['Sequence File Name', 'Format']
    TARGETLABELS = ["(>Bin) or Part Name", "Direction", "Forced Assembly Strategy?",
                    "Forced Relative Overhang Position?", "Direct Synthesis Firewall?", "Extra 5' CPEC overlap bps",
                    "Extra 3' CPEC overlap bps"]
    MASTEROLIGOSLABELS = ['Oligo Name', 'Length', "Tm", "Tm (3' only)", "Sequence"]
    MASTERPLASMIDLABELS = ["Plasmid Name", "Alias", "Contents", "Length", "Sequence"]
    MASTERDIRECTLABELS = ["Direct Synthesis Name", "Alias", "Contents", "length", "Sequence"]

    # ZIP = None
    # EUGENE = None
    # PARAMS = None

    def __init__(self, assembly):
        super(J5Assembly, self).__init__(assembly.contigs, assembly.contig_container, assembly.primer_continer)
        self.params = None
        self.proxy_home = 'https://j5.jbei.org/bin/j5_xml_rpc.pl'
        self.proxy = xmlrpclib.ServerProxy('https://j5.jbei.org/bin/j5_xml_rpc.pl')
        self.session = None
        self.session_id = None

    def to_csv(self, labels, rows):
        csvrows = []
        for row in [labels] + rows:
            csvrows.append(','.join([str(x) for x in row]))
        return '\n'.join(csvrows)

    @staticmethod
    def encode64(filestring):
        return base64.encodestring(filestring)

    @staticmethod
    def decode64(string):
        return base64.decodestring(string)

    def all(self):
        self.get_parts()
        self.get_target()
        self.get_sequences()
        self.get_plasmids()
        self.get_primers()
        self.get_direct()
        self.get_eugene()
        self.get_parameters()

    def decode_all_to(self, destination):
        self.all()
        names = [
            'parts',
            'target',
            'sequences',
            'plasmids',
            'direct',
            'eugene',
            'parameters',
        ]
        print self.__dict__.keys()
        for n in names:
            with open(os.path.join(destination, n + '.csv'), 'w') as handle:
                handle.write(J5Assembly.decode64(self.__dict__[n]))
        with open(os.path.join(destination, 'sequences.zip'), 'w') as zip:
            zip.write(J5Assembly.decode64(self.encoded_zip))

    def get_parameters(self):
        with open('j5_parameters.csv') as params:
            self.parameters = J5Assembly.encode64(params.read())

    def get_target(self):
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
        self.target = J5Assembly.encode64(self.to_csv(J5Assembly.TARGETLABELS, rows))

    def get_parts(self):
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
        self.parts = J5Assembly.encode64(csv)

    def get_direct(self):
        self.direct = J5Assembly.encode64(self.to_csv(J5Assembly.MASTERDIRECTLABELS, []))

    def get_plasmids(self):
        self.plasmids = J5Assembly.encode64(self.to_csv(J5Assembly.MASTERPLASMIDLABELS, []))

    def get_primers(self):
        # ['Oligo Name', 'Length', "Tm", "Tm (3' only)", "Sequence"]
        rows = []
        for p in self.primer_continer.contigs:
            row = [
                p.subject_acc,
                len(p.subject_seq),
                60,
                60,
                p.subject_seq
            ]
        rows.append(row)
        self.primers = J5Assembly.encode64(self.to_csv(J5Assembly.MASTEROLIGOSLABELS, rows))

    # TODO: add direct fragments to directmasterlist
    # TODO: if is direct, do something with j5 assembly parameters?
    def get_sequences(self):
        print 'SEQUENCE LIST'
        sequences = []
        for c in self.contigs:
            print '\tSEQ: {}  @  {}'.format(os.path.basename(c.filename), c.seqrecord.id)
            sequences.append((c.filename, c.seqrecord.id,))
        rows = list(set(sequences))

        # Save sequence_list

        csvrows = []
        for s in sequences:
            filename = os.path.basename(s[0])
            format = determine_format(s[0]).capitalize()
            csvrows.append((filename, format))
        csv = self.to_csv(J5Assembly.SEQLISTLABELS, csvrows)
        self.sequences = J5Assembly.encode64(csv)

        # TODO: move zip stuff to seqio
        # Save zipped sequence file
        temp = tempfile.TemporaryFile()  # open('/Users/Justin/Desktop/temporary.txt', 'w')#tempfile.TemporaryFile()
        with zipfile.ZipFile(temp, 'w') as zf:

            for s in sequences:
                filename = s[0]
                print '\tZipping {} as {}'.format(s[0], os.path.basename(filename))
                with open(filename, 'rU') as handle:
                    zf.writestr(os.path.join('ZippedTemplates', os.path.basename(filename)), handle.read())
        temp.seek(0)
        self.encoded_zip = J5Assembly.encode64(temp.read())
        temp.close()

    def get_eugene(self):
        self.eugene = J5Assembly.encode64('')

    def to_dict(self):
        d = {
            'zipped_sequences': self.encoded_zip,
            'master_oligos': self.primers,
            'master_plasmids': self.plasmids,
            'master_direct_syntheses': self.direct,
            'parts_list': self.parts,
            'j5_parameters': self.parameters,
            'sequences_list': self.sequences,
            'eugene_rules': self.eugene,
            'target_part_order_list': self.target,
        }
        params = {}
        for k in d:
            params['encoded_{}_file'.format(k)] = d[k]
            params['reuse_{}_file'.format(k)] = d[k]
        return params

    def login(self, username, password):
        self.session = self.proxy.CreateNewSessionId({'username': username, 'password': password})
        self.session_id = self.session['j5_session_id']
        return self.session_id

    def submit(self, username, password):
        self.all()
        self.login(username, password)
        return self.run_assembly()

    def run_assembly(self, assembly_method='SLIC/Gibson/CPEC'):
        d = self.to_dict()
        print d
        d['j5_session_id'] = self.session_id
        d['assembly_method'] = assembly_method
        print self.to_dict().keys()
        return self.proxy.DesignAssembly(d)

        # params['encoded_{}_file'.format(filetype)] = self.encode(imply_file('*{}*'.format(filetype)))
        #             params['reuse_{}_file'.format(filetype)] = False
