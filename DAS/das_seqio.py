'''
Project: jdna
File: das_seqio
Author: Justin
Date: 2/6/17

Description: 

'''

import shlex
import subprocess
from Bio import SeqIO
import os
import json
import coral
from glob import glob


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
    :param out: output path
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

    return out, sequences, metadata


def gb_to_fsa(input_path, output_path):
    sequences = open_sequence(input_path)
    with open(output_path, 'w') as handle:
        SeqIO.write(sequences, handle, 'fasta')
    return output_path

def seq_is_circular(path):
    d = coral.seqio.read_dna(path)
    return d.circular
