'''
Project: jdna
File: das_cost_model
Author: Justin
Date: 2/7/17

Description: 

'''
from contig import *

class CostModel:
    MIN_PCR_SIZE = 250.
    MAX_PCR_SIZE = 10000.
    PRIMER_COST = 15.
    PCR_COST = 14.
    FIVEPRIME_EXT_REACH = 20. # reachability for a extended primer
    SYNTHESIS_THRESHOLD = 60. # min bp for synthesis
    SYNTHESIS_COST = 0.11 # per bp
    SYNTHESIS_MIN_COST = 89. # per synthesis

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
        if right.end_label == 'new primer':
            r += CostModel.FIVEPRIME_EXT_REACH
        if left.start_label == 'new primer':
            r += CostModel.FIVEPRIME_EXT_REACH
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
        return query_length + CostModel.gap(right, left)

    @staticmethod
    def get_gaps(contig_path, circular=True):
        '''
        Calculates all non-circular gaps for a contig_path
        :param contig_path:
        :param circular:
        :return:
        '''
        pairs = zip(contig_path[:-1], contig_path[1:])
        gaps = []
        for l, r in pairs:
            if CostModel.assembly_condition(l, r):
                gaps.append(CostModel.gap(l, r))
            else:
                gaps.append(float("Inf"))
                return gaps
        return gaps

    @staticmethod
    def get_gap_cost(gaps):
        '''
        A monotonically increasing gap cost function.
        :param gaps:
        :return:
        '''
        non_zero_gaps = filter(lambda x: x > 0, gaps)
        return sum(non_zero_gaps)*CostModel.SYNTHESIS_COST

    @staticmethod
    def get_span_gap(contig_path, query_length):
        '''
        Get the span gap for the assembly.
        :param contig_path:
        :param query_length:
        :return:
        '''
        query_span = contig_path[-1].q_end - contig_path[0].q_start
        return query_length - query_span

    @staticmethod
    def get_span_cost(contig_path, query_length):
        return CostModel.get_span_gap(contig_path, query_length) * CostModel.SYNTHESIS_COST

    @staticmethod
    def pcr_cost(contig):
        '''
        Cost of a pcr
        :param contig:
        :return:
        '''
        if not CostModel.MIN_PCR_SIZE < contig.q_end - contig.q_start < CostModel.MAX_PCR_SIZE:
            return float("Inf")
        if not contig.circular and contig.parent == 'query':
            # direct_synthesis
            return 0
        cost = CostModel.PCR_COST
        if contig.circular:
            cost += 2 * CostModel.PRIMER_COST
        else:
            if not contig.start_label == 'primer':
                cost += CostModel.PRIMER_COST
            if not contig.end_label == 'primer':
                cost += CostModel.PRIMER_COST
        # print contig.circular, contig.parent, contig.start_label, contig.end_label, cost
        return cost

    @staticmethod
    def assembly_costs(contig_path, query_length, circular=True):
        '''
        Calculate the assembly costs for the given assembly.

        :param contig_path:
        :param query_length:
        :param circular:
        :return:
        '''
        if len(contig_path) == 0:
            return float("Inf")

        pairs = zip(contig_path[:-1], contig_path[1:])

        query_span = contig_path[-1].q_end - contig_path[0].q_start
        span_cost = query_length - query_span
        if span_cost < 0:
            span_cost = float("Inf")
        gaps = CostModel.get_gaps(contig_path)
        gap_cost = CostModel.get_gap_cost(gaps)
        span = CostModel.get_span_gap(contig_path, query_length)
        span_cost = span * CostModel.SYNTHESIS_COST
        num_synthesized_fragments = len(filter(lambda x: x > CostModel.SYNTHESIS_THRESHOLD, gaps))
        if circular:
            if span > CostModel.SYNTHESIS_THRESHOLD:
                num_synthesized_fragments += 1
            # non valid assembly

        fragment_cost = sum([CostModel.pcr_cost(c) for c in contig_path])
        return gap_cost, span_cost, num_synthesized_fragments, fragment_cost

    @staticmethod
    def assembly_probability(num_frags):
        '''
        The probability of a successful GIBSON/SLIC assembly given
        a number of fragments.
        :param num_frags:
        :return:
        '''
        p = 1.0 - (1.0 / 8.0) * num_frags
        if p <= 0:
            p = 0.0001
        return p

    @staticmethod
    def total_cost(contig_path, query_length, circular=True):
        '''
        Calculates total cost in dollars for a given assembly.
        :param contig_path:
        :param query_length:
        :param circular:
        :return:
        '''
        gap_cost, span_cost, new_frags, fragment_cost = CostModel.assembly_costs(contig_path, query_length, circular=circular)
        probability = CostModel.assembly_probability(new_frags + len(contig_path))
        return 1 / probability * (gap_cost + span_cost + fragment_cost + CostModel.SYNTHESIS_MIN_COST * new_frags)

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