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


def pcr_products_of_contig(contig, primers, ignore_direction=True, include_contig=True):
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
    return contig.divide_contig(start_positions, end_positions, include_contig=include_contig, contig_type='product')





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
        # assembly = Assembly(contigs=new_alignments)
        # a = assembly.get_all_assemblies()
        # print len(a), len(new_alignments)
        all_alignments += new_alignments
    contigs += all_alignments
    a = Assembly(contigs=contigs)
    return a.get_all_assemblies()

    # import numpy as np
    # contig_container = Assembly(contigs=contigs)
    #
    # def func(x):
    #     arr = [np.round(X) for X in x]
    #     arr = list(np.nonzero(arr)[0])
    #     d = {i: c for i, c in enumerate(contig_container.contigs)}
    #     contigs = [d[x] for x in arr]
    #     contigs = sorted(contigs, key=lambda x: x.q_start)
    #     cost = contig_container.compute_assembly_cost(contigs)
    #     print cost
    #     return cost
    # # contig_graph = create_alignment_graph(new_alignments)
    # # contig_graph['root'] = [x['align_id'] for x in new_alignments]
    # # paths = dfs_iter(contig_graph, 'root', new_alignments)
    # # a['assemblies'] = paths
    # from pyswarm import pso
    # lb = [0] * len(contig_container.contigs)
    # ub = [1] * len(contig_container.contigs)
    # xopt, fopt = pso(func, lb, ub, swarmsize=100, maxiter=100)
    # arr = [np.round(X) for X in xopt]
    # arr = list(np.nonzero(arr)[0])
    # d = {i: c for i, c in enumerate(contig_container.contigs)}
    # contigs = [d[x] for x in arr]
    # print [(c.q_start, c.q_end) for c in contigs]

        # TODO: suggestions for splitting large pcr products
        # TODO: trim graph if path spans the query length
        # TODO: cost calculation for distance betweeen last contig and query_end?
