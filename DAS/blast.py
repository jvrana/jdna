from glob import glob
import subprocess
import shlex
from Bio import SeqIO
import os

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

def concat_seqs(idir, odir, out_name):
    sequences = []

    formats = {"genbank": ["gb"], "fasta": ["fsa", "fasta"]}

    for fmt in formats:
        filenames = []
        for suffix in formats[fmt]:
            paths = os.path.join(idir, '*.{}'.format(suffix))
            filenames += glob(paths)

        for filename in filenames:
            with open(filename, 'rU') as input_handle:
                sequences += SeqIO.parse(input_handle, fmt)

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
    sequences = []
    with open(input, 'rU') as input_handle, open(output, 'w') as output_handle:
        sequences += SeqIO.parse(input_handle, 'genbank')
        SeqIO.write(sequences, output_handle, 'fasta')
    return output

def runblast(db, query, out, outfmt=10):
    "blastn -db db.fsa -query 9320_\ pLSL-[ADH1]-p35S\:mTURQ\:t35S-[ADH1]-R.fsa -out results.out"

    prefix, suffix = query.split('.')
    if suffix == 'gb':
        query = gb_to_fsa(query, prefix + '.fsa')

    cmd_str = "blastn -db {db} -query {query} -out {out} -outfmt {outfmt}"
    run_cmd(cmd_str, db=db, query=query, out=out, outfmt=outfmt)
    return out

c = concat_seqs(locations['templates'], locations['database'], 'db')
db = makedb(c, 'db', os.path.join(locations['database'], 'db'))



runblast(db,
         os.path.join(locations['designs'], 'pmodkan-ho-pact1-z4-er-vpr.gb'),
         os.path.join(locations['database'], 'results.out')
         )
# def blastn_cmd(db='', query='', out=''):
#     cmd_str = "blastn -db {database} -query {query} -out {out}".format(database=db, query=query, out=out)
#     args = shlex.strip(cmd_str)
#     return args

