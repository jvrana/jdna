'''
Project: jdna
File: das_assembler
Author: Justin
Date: 2/6/17

Description: 

'''

import itertools
from collections import defaultdict
from das_blast import *
from copy import deepcopy
from contig import *


def dfs_iter(graph, root, alignments):
    paths = []
    stack = [[root]]
    best_costs = [float("Inf")] * 10  # five best costs
    while stack:
        best_costs.sort()
        path = stack.pop()
        cost = compute_cost(path, alignments)
        # Trimming conditions
        if cost >= best_costs[-1] and not cost == float("Inf"):
            # cost can only get worst, so trim this search
            continue

        # if path > query_length

        n = path[-1]
        if n in graph:
            for node in graph[n]:
                new_path = path[:] + [node]
                stack.append(new_path)
        else:
            if cost < best_costs[-1]:
                best_costs[-1] = cost
                print cost
            paths.append(path)
    for path in paths:
        path.remove('root')
    return paths


# TODO: add pcr length cost
# TODO: add cost calculation for each fragment
# TODO: must be monotonically increasing
# TODO: add cost if primer doesn't exist





# Todo: Double all dnas to account for cyclic plasmids (if cyclic)
# TODO: if linearized, can use directly
# TODO: find all primer positions
# TODO: account for reversed alignments



def get_primers_within_bounds(contig, primers, minimum_primer_anneal=15):
    def primer_within_bounds(primer):
        minimum_primer_anneal = 15  # TODO: move MIN_PRIMER to design parameters
        if primer.subject_strand == 'plus':
            return contig.q_start + minimum_primer_anneal < primer.q_end < contig.q_end
        elif primer.subject_strand == 'minus':
            return contig.q_start < primer.q_start < contig.q_end - minimum_primer_anneal
        return False

    return filter(lambda x: primer_within_bounds(x), primers)


def pcr_products_of_contig(contig, primers, ignore_direction=True):
    '''
    Produces all possible fragments from primer positions
    :param contig:
    :param primers:
    :param ignore_direction:
    :return:
    '''

    new_alignments = []
    primers = get_primers_within_bounds(contig, primers)
    end_positions = []
    start_positions = []
    for primer in primers:
        direction = primer.subject_strand
        if direction == 'plus' or ignore_direction:
            start_positions.append(primer.q_start)
        if direction == 'minus' or ignore_direction:
            end_positions.append(primer.q_end)
    return contig.divide_contig(start_positions, end_positions, includecontig=True, contig_type='product')


def compute_cost(contig_path):
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


def assembly_condition(left, right):
    r_5prime_threshold = 0
    r_3prime_threshold = 60
    l_pos = left.q_end
    r_pos = right.q_start
    # return r_pos == l_pos + 1 # if its consecutive
    return r_pos > l_pos - r_3prime_threshold and \
           right.q_end > left.q_end


def assembly_graph(contig_container, sort=True):
    assembly_graph = {}
    pairs = list(itertools.permutations(contig_container.contigs, 2))
    print '{} pairs to search'.format(len(pairs))
    for l, r in pairs:

        if l == r:
            continue
        if assembly_condition(l, r):
            if l.contig_id not in assembly_graph:
                assembly_graph[l.contig_id] = []
            assembly_graph[l.contig_id].append(r.contig_id)

    if sort:
        for k in assembly_graph.keys():
            a_array = assembly_graph[k][:]

            assembly_graph[k] = sorted(a_array, key=lambda x: compute_cost([contig_container.get_contig(k),contig_container.get_contig(x)]))
    return assembly_graph
    # Find Paths


# TODO: generate additional pcr products that could be homologous to existing primer
# TODO: compute alignment_graph for each 'contig'
def contig_assembly(contigs, primers):
    '''
    Calculates paths for contig assemblies
    :param alignment_results:
    :param primer_results:
    :return:
    '''
    all_alignments = []
    for contig in contigs:
        new_alignments = pcr_products_of_contig(contig, primers)
        all_alignments += new_alignments
    contigs += all_alignments
        # contig_graph = create_alignment_graph(new_alignments)
        # contig_graph['root'] = [x['align_id'] for x in new_alignments]
        # paths = dfs_iter(contig_graph, 'root', new_alignments)
        # a['assemblies'] = paths

        # TODO: suggestions for splitting large pcr products
        # TODO: trim graph if path spans the query length
        # TODO: cost calculation for distance betweeen last contig and query_end?
