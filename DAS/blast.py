from glob import glob
import subprocess
import shlex
from Bio import SeqIO
import os
import re
import json

locations = dict(
    database="database",
    designs="designs",
    templates="templates"
)

"makeblastdb -in db.fsa -dbtype nucl -title name -out"

"blastn -db db.fsa -query 9320_\ pLSL-[ADH1]-p35S\:mTURQ\:t35S-[ADH1]-R.fsa -out results.out"


def run_cmd(cmd_str, **kwargs):
    args = shlex.split(cmd_str.format(**kwargs))
    output = subprocess.Popen(args)

def open_sequence(filename, fullname=True, **fmt):
    formats = {"gb": "genbank", "fsa": "fasta", "fasta": "fasta"}
    formats.update(fmt)
    prefix, suffix = os.path.basename(filename).split('.')
    seqs = []
    with open(filename, 'rU') as output_handle:
        s = list(SeqIO.parse(output_handle, formats[suffix]))
        # Rename accession id to full name
        if len(s) == 1 and fullname:
            s[0].id = prefix
        seqs += s
    return seqs

def concat_seqs(idir, odir, out_name):
    sequences = []

    filenames = glob(os.path.join(idir, '*.*'))
    for filename in filenames:
        print 'opening', filename
        sequences += open_sequence(filename)
    output_path = os.path.join(odir, out_name + '.fsa')
    with open(output_path, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")
    return output_path

def makedb(input, title, output):
    sequences = []
    for filename in glob(os.path.join(locations['templates'], '*.gb')):
        with open(filename, "rU") as input_handle:
            sequences += SeqIO.parse(input_handle, "genbank")

    with open("allsequences.fasta", "w") as output_handle:
         SeqIO.write(sequences, output_handle, "fasta")


    cmd_str = "makeblastdb -in {input} -dbtype nucl -title {title} -out {output}"
    run_cmd(cmd_str, input=input, title=title, output=output)
    return output

def gb_to_fsa(input, output):
    sequences = open_sequence(input)
    with open(output, 'w') as output_handle:
        SeqIO.write(sequences, output_handle, 'fasta')
    return output

def runblast(db, query, out, output_format=7, output_views="qacc sacc score evalue bitscore\
 length nident gapopen gaps qlen qstart qend slen sstart send sstrand qseq sseq"):

    prefix, suffix = query.split('.')
    if suffix == 'gb':
        query = gb_to_fsa(query, prefix + '.fsa')
    outfmt = '\"{} {}\"'.format(output_format, output_views)
    cmd_str = "blastn -db {db} -query {query} -out {out} -outfmt {outfmt}"
    run_cmd(cmd_str, db=db, query=query, out=out, outfmt=outfmt)

    results = ''
    with open(out, 'rU') as output_handle:
        results = output_handle.read()

    return results
# Query: pMODKan-HO-pACT1-out

c = concat_seqs(locations['templates'], locations['database'], 'db')
db = makedb(c, 'db', os.path.join(locations['database'], 'db'))
# query acc., subject acc., % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score


def parse_results(results):
    g = re.search('# (?P<blast_ver>.+)\n# Query: (?P<query>.+)\n# Database: (?P<database>.+)\n# Fields: (?P<fields>.+)', results)
    print g
    print results
    meta = g.groupdict()

    matches = re.findall('\n([^#].*)', results)
    meta['fields'] = re.split('\s*,\s*', meta['fields'])
    # clean up fields
    for i, f in enumerate(meta['fields']):
        meta['fields'][i] = f.replace('.', '')\
            .replace(' ', '_')\
            .replace('%', 'perc')
    alignments = []
    for m in matches:
        alignments.append(
            dict(zip(meta['fields'], m.split('\t')))
        )

    return {'meta': meta, 'alignments': alignments}

# Todo: Double all dnas to account for cyclic plasmids (if cyclic)
# TODO: if linearized, can use directly
# TODO: find all primer positions

r = runblast(db,
         os.path.join(locations['designs'], 'pmodkan-ho-pact1-z4-er-vpr_doubled_temp.gb'),
         os.path.join(locations['database'], 'results.out'))
p = parse_results(r)
with open('alignment_viewer/data.json', 'w') as output_handle:
    json.dump(p, output_handle)
# def blastn_cmd(db='', query='', out=''):
#     cmd_str = "blastn -db {database} -query {query} -out {out}".format(database=db, query=query, out=out)
#     args = shlex.strip(cmd_str)
#     return args

