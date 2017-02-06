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

def run_cmd(cmd_str, **kwargs):
    """
    Runs a subprocess command with kwargs arguments
    :param cmd_str: formattable command line string (e.g. 'ls {directory_location}')
    :param kwargs: dictionary of arguments for command
    :return: None
    :param cmd_str:
    :param kwargs:
    :return:
    """
    cmd_str = cmd_str.format(**kwargs)
    args = shlex.split(cmd_str)
    print "CMD: {}".format(cmd_str)
    output = subprocess.Popen(args)
    output.wait()


def format_decorator(f):
    """
    Implies a SeqIO format based on filename suffix
    :param f:
    :return: wrapped function
    """

    def wrapped(*args, **kwargs):
        args = list(args)
        prefix, suffix = args[0].split('.')
        formats = {"gb": "genbank", "fsa": "fasta", "fasta": "fasta"}
        formats.update(kwargs)
        kwargs['format'] = formats[suffix]
        return f(*args, **kwargs)

    return wrapped

def assign_id(alignment):
    global ALIGN_ID
    ALIGN_ID += 1
    alignment['align_id'] = ALIGN_ID

@format_decorator
def open_sequence(filename, format=None, **fmt):
    prefix, suffix = os.path.basename(filename).split('.')
    seqs = []
    with open(filename, 'rU') as handle:
        s = list(SeqIO.parse(handle, format))
        if len(s) == 1:
            s[0].id = prefix
        seqs += s
    return seqs


@format_decorator
def save_sequence(filename, sequences, format=None, **fmt):
    with open(filename, 'w') as handle:
        SeqIO.write(sequences, handle, format)
    return filename


def concat_seqs(idir, out, savemeta=False):
    """
    Concatenates a directory of sequences into a single fasta file
    :param idir: input directory path
    :param out: output path (no suffix)
    :return: full output path
    """
    sequences = []
    metadata = {}
    filenames = glob(os.path.join(idir, '*.*'))
    for filename in filenames:
        seqs = open_sequence(filename)
        sequences += seqs

        # TODO: this is really hacky, recode this
        for s in seqs:
            metadata[s.id] = {'filename': filename, 'circular': False}
        if len(seqs) == 1:
            c = coral.seqio.read_dna(filename)
            metadata[seqs[0].id]['circular'] = c.circular

    with open(out, "w") as handle:
        SeqIO.write(sequences, handle, "fasta")

    if savemeta:
        with open(out.split('.')[0] + '.json', 'w') as handle:
            json.dump(metadata, handle)

    return out, sequences


def makedb(fasta, title, database_output):
    """
    Makes a blast database from a fasta file
    :param fasta: path to input fasta file
    :param title: title to name database (e.g. "db")
    :param database_output: path to output directory (e.g. "/Users/databases/db")
    :return: path to output directory
    """
    cmd_str = "makeblastdb -in {input} -dbtype nucl -title {title} -out {output}"
    run_cmd(cmd_str, input=fasta, title=title, output=database_output)

    return database_output


def makedbfromdir(idir, out):
    fasta, seqs = concat_seqs(idir, out + '.fsa', savemeta=True)
    makedb(fasta, os.path.basename(out), out)
    return out


def gb_to_fsa(input_path, output_path):
    sequences = open_sequence(input_path)
    with open(output_path, 'w') as handle:
        SeqIO.write(sequences, handle, 'fasta')
    return output_path


def runblast(db, query, out, output_format=7, output_views="qacc sacc score evalue bitscore\
 length nident gapopen gaps qlen qstart qend slen sstart send sstrand qseq sseq", **additional_params):
    prefix, suffix = query.split('.')
    if suffix == 'gb':
        query = gb_to_fsa(query, prefix + '.fsa')
    outfmt = '\"{} {}\"'.format(output_format, output_views)
    cmd_str = "blastn -db {db} -query {query} -out {out} -outfmt {outfmt}"
    for key in additional_params:
        cmd_str += ' -{} {}'.format(key, additional_params[key])
    run_cmd(cmd_str, db=db, query=query, out=out, outfmt=outfmt)

    with open(out, 'rU') as handle:
        return handle.read()


def parse_results(results, additional_metadata=None):
    g = re.search(
        '# (?P<blast_ver>.+)\n# Query: (?P<query>.+)\\n# Database: (?P<database>.+)\n# Fields: (?P<fields>.+)', results)

    meta = g.groupdict()

    matches = re.findall('\n([^#].*)', results)
    meta['fields'] = re.split('\s*,\s*', meta['fields'])
    # clean up fields
    for i, f in enumerate(meta['fields']):
        meta['fields'][i] = f.replace('.', '') \
            .replace(' ', '_') \
            .replace('%', 'perc')
    alignments = []
    for m in matches:
        v = m.split('\t')

        # infer datatypes
        for i, val in enumerate(v):
            try:
                v[i] = float(val)
            except ValueError as e:
                pass
            try:
                v[i] = int(val)
            except ValueError as e:
                pass

        align_dict = dict(zip(meta['fields'], v))
        assign_id(align_dict)
        align_dict['align_type'] = 'alignment'
        alignments.append(
            align_dict
        )
        if additional_metadata is not None:
            if align_dict['subject_acc'] in additional_metadata:
                align_dict.update(additional_metadata[align_dict['subject_acc']])
    return {'meta': meta, 'alignments': sorted(alignments, key=lambda x: x['q_start'])}


def fill_in_gaps(alignments):
    pairs = itertools.permutations(alignments, 2)
    gaps = []
    for l, r in pairs:
        l_end = l['q_end']
        r_start = r['q_start']
        if l_end < r_start:
            start = l_end + 1
            end = r_start - 1
            if start <= end:
                g = {'q_start': start, 'q_end': end, 'align_type': 'gap'}
                assign_id(g)
                gaps.append(g)
    print '{} gaps'.format(len(gaps))
    alignments += sorted(gaps, key=lambda x: x['q_start'])


def alignment_condition(left, right):
    r_5prime_threshold = 0
    r_3prime_threshold = 60
    l_pos = left['q_end']
    r_pos = right['q_start']
    # return r_pos == l_pos + 1 # if its consecutive
    return r_pos > l_pos - r_3prime_threshold and \
            right['q_end'] > left['q_end']



def create_alignment_graph(alignments):
    alignment_graph = defaultdict(list)
    # alignments = p['alignments']
    pairs = list(itertools.permutations(alignments, 2))
    print '{} pairs to search'.format(len(pairs))
    for l, r in pairs:
        if l == r:
            continue
        if alignment_condition(l, r):
            alignment_graph[l['align_id']].append(r['align_id'])

    return dict(alignment_graph)
    # Find Paths


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

def make_primer_alignments():
    seq = open_sequence('designs/pmodkan-ho-pact1-z4-er-vpr.gb')[0]
    import numpy as np
    min_size = 15
    max_size = 60
    num_primers = 10
    rand_pos = np.random.randint(0, len(seq)-max_size, size=num_primers)
    rand_size = np.random.randint(min_size, max_size, size=num_primers)
    primers = []
    for pos, size in zip(rand_pos, rand_size):
        primer = seq[pos:pos+size]
        primer.id = primer.id + '{}-{}'.format(pos, pos+size)

        rc = np.random.choice(['rc', ''])
        if rc == 'rc':
            primer = primer.reverse_complement(id=primer.id+'_rc')
        primers.append(primer)
    save_sequence("primers/primers.fasta", primers)

    # TODO: add new alignments to p['alignments'], these cost much more (fragment_construction_cost)


    database = makedbfromdir('primers', os.path.join(locations['database'], 'primerdb'))

    # Run blast
    r = runblast(database,
                 'designs/newseq.gb',
                 os.path.join(locations['database'], 'primer_results.out'), word_size=15, evalue=1000)

    p = parse_results(r)
    with open('alignment_viewer/primer_data.json', 'w') as output_handle:
        json.dump(p, output_handle)
    return p

def run_alignment():
    print 'Running alignment against templates'
    # Make database
    database = makedbfromdir(locations['templates'], os.path.join(locations['database'], 'db'))

    # TODO: only if circular, then pseudo-circular
    seqs = open_sequence(os.path.join(locations['designs'], 'pmodkan-ho-pact1-z4-er-vpr.gb'))
    seq = seqs[0]
    save_sequence('designs/newseq.gb', seq + seq)

    # Run blast
    r = runblast(database,
                 'designs/newseq.gb',
                 os.path.join(locations['database'], 'results.out'), evalue=10)

    # fuse circular alignments

    # open meta data
    meta = {}
    with open('database/db.json', 'rU') as handle:
        meta = json.load(handle)

    # parse results
    p = parse_results(r, additional_metadata=meta)






    # clean up alignments
    fuse_circular_fragments(p['alignments'])

    primer_results = make_primer_alignments()
    # TODO: filter out primers with imperfect binding
    # TODO: filter out alignments with imperfect binding

    def primer_within_bounds(primer, alignment):
        MIN_PRIMER = 15 # TODO: move MIN_PRIMER to design parameters
        if primer['subject_strand'] == 'plus':
            # print alignment['q_start'] + MIN_PRIMER, primer['q_end'], alignment['q_end'], alignment['q_start'] + MIN_PRIMER < primer['q_end'] < alignment['q_end']
            return alignment['q_start'] + MIN_PRIMER < primer['q_end'] < alignment['q_end']
        elif primer['subject_strand'] == 'minus':
            return alignment['q_start'] < primer['q_start'] < alignment['q_end'] - MIN_PRIMER
        return False

    new_alignments = []
    for a in p['alignments']:
        primers = filter(lambda x: primer_within_bounds(x, a), primer_results['alignments'])
        end_positions = [a['q_end']]
        start_positions = [a['q_start']]
        for primer in primers:
            direction = primer['subject_strand']
            if direction == 'plus':
                start_positions.append(primer['q_start'])
            elif direction == 'minus':
                end_positions.append(primer['q_end'])
        new_positions = list(itertools.product(start_positions, end_positions))
        for s, e in new_positions:
            if s >= e:
                continue
            new_a = copy(a)
            new_a['q_start'] = s
            new_a['q_end'] = e
            if 's_start' in new_a:
                new_a['s_start'] = a['s_start'] + (s - a['q_start'])
                new_a['s_end'] = a['s_end'] - (a['q_end'] - e)
            new_a['alignment_length'] = e - s
            assign_id(new_a)
            new_a['align_type'] = 'product'
            new_alignments.append(new_a)

    print 'sorting primer alignments'
    p['alignments'] += sorted(new_alignments, key=lambda x: x['q_start'])

    print 'filling in gaps'
    # fill_in_gaps(p['alignments'])
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

# compute cost of path as graph is built

# collect all contigs, calculate cost of gap

print 'creating alignment_graph'
alignment_graph = create_alignment_graph(p['alignments'])


d = {a['align_id']: a for a in p['alignments']}



#
def compute_cost(path, alignments):
    d = {a['align_id']: a for a in alignments}
    if path[0] not in d:
        return float("Inf")
    first = d[path[0]]

    cost = first['q_start']

    pairs = zip(path[:-1], path[1:])
    L, R = -1, -1
    for l, r in pairs:
        L = d[l]
        R = d[r]
        if first['q_start'] - R['q_end'] > alignments[0]['query_length']/2.0:
            continue
        cost += R['q_start'] - L['q_end']
    if R == -1:
        return float("Inf")
    # cost += alignments[0]['query_length'] - R['q_end']
    return cost


def dfs_iter(graph, root):
    paths = []
    stack = [[root]]
    best_costs = [float("Inf")]*5 # five best costs
    while stack:
        best_costs.sort()
        path = stack.pop()
        cost = compute_cost(path, p['alignments'])
        if cost > best_costs[-1]:
            # cost can only get worst, so trim this search
            continue
        n = path[-1]
        if n in graph:
            for node in graph[n]:
                new_path = path[:] + [node]
                stack.append(new_path)
        else:
            if cost < best_costs[-1]:
                best_costs[-1] = cost
            paths.append(path)
    return paths

# alignment_graph['root'] = [a['align_id'] for a in p['alignments']]
paths = []
for i, a in enumerate(p['alignments']):
    print i, len(p['alignments'])
    paths += dfs_iter(alignment_graph, a['align_id'])

paths = sorted(paths, key=lambda x: compute_cost(x, p['alignments']))
print paths
print len(paths)
print [compute_cost(path, p['alignments']) for path in paths]
best = paths[:5]
d = {a['align_id']: a for a in p['alignments']}
for path in best:
    print [(d[id]['q_start'], d[id]['q_end']) for id in path], compute_cost(path, p['alignments'])
#
d = {a['align_id']: a for a in p['alignments']}
ids = [x['align_id'] for x in p['alignments']]
path = paths[0]
print
assert len(ids) == len(set(ids))

# Compute costs

# run j5
# TODO: make sure reverse complements are handled correctly

# compute costs

# display alignments

# # Save results


