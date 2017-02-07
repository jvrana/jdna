from glob import glob
import subprocess
import shlex
from Bio import SeqIO
import os
import re
import json
from collections import defaultdict
import coral
import itertools
from copy import copy

locations = dict(
    database="database",
    designs="designs",
    templates="templates"
)

"makeblastdb -in db.fsa -dbtype nucl -title name -out"

"blastn -db db.fsa -query 9320_\ pLSL-[ADH1]-p35S\:mTURQ\:t35S-[ADH1]-R.fsa -out results.out"


ALIGN_ID = 0






def _recursive_dfs(graph, start, path=[]):
    '''recursive depth first search from start'''
    path = path + [start]
    if start in graph:
        for node in graph[start]:
                if not node in path:
                    path = _recursive_dfs(graph, node, path)
    return path


def recursive_dfs(graph):
    paths = []
    for start in graph:
        paths.append(_recursive_dfs(graph, start, []))
    return paths


def fuse_circular_fragments(alignments):
    ga = defaultdict(list)

    def fuse_condition(l, r):
        return l['subject_acc'] == r['subject_acc'] and \
               l['q_end'] + 1 == r['q_start'] and \
               l['subject_length'] == l['s_end'] and \
               l['circular'] and r['circular']

    def fuse(l, r):
        l['q_end'] = r['q_end']
        l['s_end'] = r['s_end']
        keys_to_sum = 'alignment_length, identical, gap_opens, gaps, identical, score'.split(', ')
        for k in keys_to_sum:
            l[k] += r[k]

    pairs = itertools.permutations(alignments, 2)
    for l, r in pairs:
        if fuse_condition(l, r):
            fuse(l, r)
            alignments.remove(r)


# Todo: Double all dnas to account for cyclic plasmids (if cyclic)
# TODO: if linearized, can use directly
# TODO: find all primer positions
# TODO: account for reversed alignments

def generate_random_primers(seq, out):
    import numpy as np
    min_size = 15
    max_size = 60
    num_primers = 5
    rand_pos = np.random.randint(0, len(seq) - max_size, size=num_primers)
    rand_size = np.random.randint(min_size, max_size, size=num_primers)
    primers = []
    for pos, size in zip(rand_pos, rand_size):
        primer = seq[pos:pos + size]
        primer.id = primer.id + '{}-{}'.format(pos, pos + size)

        rc = np.random.choice(['rc', ''])
        if rc == 'rc':
            primer = primer.reverse_complement(id=primer.id + '_rc')
        primers.append(primer)
    save_sequence(out, primers)
    return out

def make_primer_alignments(primer_dir, template_path):
    # TODO: add new alignments to p['alignments'], these cost much more (fragment_construction_cost)


    database = makedbfromdir('primers', os.path.join(locations['database'], 'primerdb'))

    # Run blast
    r = runblast(database,
                 template_path,
                 os.path.join(locations['database'], 'primer_results.out'), word_size=15, evalue=1000)

    p = parse_results(r)
    with open('alignment_viewer/primer_data.json', 'w') as output_handle:
        json.dump(p, output_handle)
    return p

def produce_fragments_from_contig(alignment, primers, ignore_direction=True):
    '''
    Produces all possible fragments from primer positions
    :param alignment:
    :param primers:
    :param ignore_direction:
    :return:
    '''
    def primer_within_bounds(primer):
        MIN_PRIMER = 15  # TODO: move MIN_PRIMER to design parameters
        if primer['subject_strand'] == 'plus':
            # print alignment['q_start'] + MIN_PRIMER, primer['q_end'], alignment['q_end'], alignment['q_start'] + MIN_PRIMER < primer['q_end'] < alignment['q_end']
            return alignment['q_start'] + MIN_PRIMER < primer['q_end'] < alignment['q_end']
        elif primer['subject_strand'] == 'minus':
            return alignment['q_start'] < primer['q_start'] < alignment['q_end'] - MIN_PRIMER
        return False

    new_alignments = []
    primers = filter(lambda x: primer_within_bounds(x), primers)
    end_positions = [alignment['q_end']]
    start_positions = [alignment['q_start']]
    for primer in primers:
        direction = primer['subject_strand']
        if direction == 'plus' or ignore_direction:
            start_positions.append(primer['q_start'])
        if direction == 'minus' or ignore_direction:
            end_positions.append(primer['q_end'])
    new_positions = list(itertools.product(start_positions, end_positions))
    for s, e in new_positions:
        if s >= e:
            continue
        new_a = copy(alignment)
        new_a['q_start'] = s
        new_a['q_end'] = e
        if 's_start' in new_a:
            new_a['s_start'] = alignment['s_start'] + (s - alignment['q_start'])
            new_a['s_end'] = alignment['s_end'] - (alignment['q_end'] - e)
        new_a['alignment_length'] = e - s
        assign_id(new_a)
        new_a['align_type'] = 'product'
        # if not (s == a['q_start'] and e == a['q_end']):
        new_alignments.append(new_a)

# TODO: generate additional pcr products that could be homologous to existing primer
# TODO: compute alignment_graph for each 'contig'
def contig_assembly(alignment_results, primer_results):
    '''
    Calculates paths for contig assemblies
    :param alignment_results:
    :param primer_results:
    :return:
    '''

    old_alignments = alignment_results['alignments']
    for a in old_alignments:
        new_alignments = produce_fragments_from_contig(a, primer_results['alignments'])

        contig_graph = create_alignment_graph(new_alignments)
        contig_graph['root'] = [x['align_id'] for x in new_alignments]
        paths = dfs_iter(contig_graph, 'root', new_alignments)
        a['assemblies'] = paths

# TODO: suggestions for splitting large pcr products
# TODO: trim graph if path spans the query length
# TODO: cost calculation for distance betweeen last contig and query_end?
def run_alignment():
    print 'Running alignment against templates'
    # Make database
    database = makedbfromdir(locations['templates'], os.path.join(locations['database'], 'db'))

    # TODO: only if circular, then pseudo-circular

    circular = True
    design_seq = os.path.join(locations['designs'], 'pmodkan-ho-pact1-z4-er-vpr.gb')
    seq = open_sequence(design_seq)[0]
    if circular:
        prefix, suffix = design_seq.split('.')
        design_seq = prefix + '_pseudocircular.' + suffix
        save_sequence(design_seq, seq + seq)

    # Run blast
    r = runblast(database,
                 design_seq,
                 os.path.join(locations['database'], 'results.out'), evalue=10)

    # fuse circular alignments

    # open meta data
    meta = {}
    with open('database/db.json', 'rU') as handle:
        meta = json.load(handle)

    # parse results
    p = parse_results(r, additional_metadata=meta)
    p['meta']['query_circular'] = circular
    p['meta']['query_length'] = len(seq)

    # clean up alignments
    fuse_circular_fragments(p['alignments'])



    generate_random_primers(seq, 'primers/primers.fasta')
    primer_results = make_primer_alignments(seq, design_seq)
    generate_pcr_products(p, primer_results)

    # TODO: filter out primers with imperfect binding
    # TODO: filter out alignments with imperfect binding

    '''
    for a in alignments:
        fwd_primer_positions = []
        rev_primer_positions = []
        def within_bounds(p):
            return if_within_bounds

        start_positions = [q_start] + forward primer starts within bounds
        end_positions = [q_end] + reverse primer starts within bounds
        itertools.product(start_positions, end_positions) - [q_start, q_end]
        align_type = 'potential pcr?'
    '''

    # p['alignments'] = sorted(p['alignments'], key=lambda x: x['q_start'])
    # Save results
    print 'saving results ({})'.format(len(p['alignments']))
    with open('alignment_viewer/data.json', 'w') as output_handle:
        json.dump(p, output_handle)

    return p


p = run_alignment()

# fill_in_gaps(p['alignments'])



# Primers

# eliminate search space for filtering by gap size

# eliminate redundent alignments

# eliminate imperfect alignments

# eliminate imperfect primer bindings

#  cost of path as graph is built

# collect all contigs, calculate cost of gap



# def create_graph()
#
# # alignment_graph['root'] = [a['align_id'] for a in p['alignments']]
# paths = []
#
# print 'creating alignment_graph'
# alignment_graph = create_alignment_graph(p['alignments'])
#
# d = {a['align_id']: a for a in p['alignments']}
#
# d = {a['align_id']: a for a in p['alignments']}
# alignment_graph['root'] = alignment_graph.keys()
#
# # exit()
# # for i, a in enumerate(p['alignments']):
# #     print i, len(p['alignments'])
# paths = dfs_iter(alignment_graph, 'root')
#
# paths = sorted(paths, key=lambda x: compute_assembly_cost(x, p['alignments']))
#
# print paths
# print len(paths)
# print [compute_assembly_cost(path, p['alignments']) for path in paths]
# best = paths[:5]
#
# for path in best:
#     print [(d[id]['q_start'], d[id]['q_end']) for id in path], compute_assembly_cost(path, p['alignments'])
# #
# d = {a['align_id']: a for a in p['alignments']}
# ids = [x['align_id'] for x in p['alignments']]
# path = paths[0]
# print
# assert len(ids) == len(set(ids))

# Compute costs

# run j5
# TODO: make sure reverse complements are handled correctly

# compute costs

# display alignments

# # Save results


