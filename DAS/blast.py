from glob import glob
import subprocess
import shlex
from Bio import SeqIO
import os
import re
import json
from collections import defaultdict
import coral

locations = dict(
    database="database",
    designs="designs",
    templates="templates"
)

"makeblastdb -in db.fsa -dbtype nucl -title name -out"

"blastn -db db.fsa -query 9320_\ pLSL-[ADH1]-p35S\:mTURQ\:t35S-[ADH1]-R.fsa -out results.out"


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


def parse_results(results):
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


        alignments.append(
            dict(zip(meta['fields'],v))
        )

    return {'meta': meta, 'alignments': alignments}


# Todo: Double all dnas to account for cyclic plasmids (if cyclic)
# TODO: if linearized, can use directly
# TODO: find all primer positions
# TODO: account for reversed alignments

# Make database
database = makedbfromdir(locations['templates'], os.path.join(locations['database'], 'db'))


# pseudo-circular
seqs = open_sequence(os.path.join(locations['designs'], 'pmodkan-ho-pact1-z4-er-vpr.gb'))
seq = seqs[0]
save_sequence('designs/newseq.gb', seq+seq)


# Run blast
r = runblast(database,
             'designs/newseq.gb',
             os.path.join(locations['database'], 'results.out'), evalue=10)
p = parse_results(r)



ga = defaultdict(list) # alignments grouped by subject_id

for a in p['alignments']:
    i = a['subject_acc']
    ga[i].append(a)

import itertools

# fuse circular alignments
meta = {}
with open('database/db.json', 'rU') as handle:
    meta = json.load(handle)

for key in ga:

    alignments = ga[key]
    if not meta[key]['circular']:
        continue
    # get only perfect alignments
    # remove duplicates (starting from back?)
    # fuse alignments
        # if subject is circular
        # if subject_end == length and q_end and q_start are consecutive
    pairs = itertools.permutations(alignments, 2)
    pairs_to_fuse = []
    for p1, p2 in pairs:
        if int(p1['s_end']) == int(p1['subject_length']):
            if int(p1['q_end']) + 1 == int(p2['q_start']):
                pairs_to_fuse.append((p1, p2))

    alignments_to_remove = []
    if len(pairs_to_fuse) == 1:
        p1, p2 = pairs_to_fuse[0]
        p1['q_end'] = p2['q_end']
        p1['s_end'] = p2['s_end']
        keys_to_sum = 'alignment_length, identical, gap_opens, gaps, identical, score'.split(', ')
        for k in keys_to_sum:
            p1[k] = p1[k] + p2[k]

        for a in p['alignments']:
            if a['s_start'] == p2['s_start']:
                alignments_to_remove.append(a)
        for a in alignments_to_remove:
            p['alignments'].remove(a)

        # TODO: remove redunent alignments?
        # remove all p1 alignments
        # remove all p2 alignments

        # print new_p
    elif len(pairs_to_fuse) > 1:
        raise Exception("There is more than one alignment pair to fuse. There must be overlapping alignments.")


# TODO: add new alignments to p['alignments'], these cost much more (fragment_construction_cost)

alignment_graph = defaultdict(list)
x1 = 60
x2 = 60

alignments = p['alignments']
for a in alignments:
    for a2 in alignments:
        if a == a2:
            continue
        x = a['q_end']
        if x - 60 <= a2['q_start'] <= x + 60:
            alignment_graph[a['q_start']].append(a2['q_start'])

# Find Paths

def recursive_dfs(graph, start, path=[]):
  '''recursive depth first search from start'''
  path=path+[start]
  for node in graph[start]:
    if not node in path:
      path=recursive_dfs(graph, node, path)
  return path

print recursive_dfs(alignment_graph,alignments[0]['q_start'])

# Compute costs

# run j5

# compute costs

# display alignments

# # Save results
with open('alignment_viewer/data.json', 'w') as output_handle:
    json.dump(p, output_handle)
