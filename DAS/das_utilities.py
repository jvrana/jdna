'''
Project: jdna
File: das_utilities
Author: Justin
Date: 2/6/17

Description: 

'''

from das_seqio import *
import numpy as np
import coral
import json



def generate_random_primers(seq, out, num_primers=5, min_size=15, max_size=60):
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

def dump_coral_to_json(path, outpath, width=700):
    seq = coral.seqio.read_dna(path)
    print path
    print seq.features
    print seq.features[0].__dict__
    print seq.__dict__.keys()

    def to_json(c, outpath):
        d = c.__dict__
        d['width'] = width
        d['height'] = d['width']
        d['length'] = len(c)
        d['features'] = [x.__dict__ for x in c.features]
        del d['top']
        del d['bottom']
        with open(outpath, 'w') as handle:
            json.dump(d, handle)

    to_json(seq, outpath)

def dna_complement(sequence):
    d1 = 'atgcn'
    d2 = 'tacgn'
    dic = dict(
        zip(
            list(d1.lower()) + list(d1.upper()),
            list(d2.lower()) + list(d2.upper())
        )
    )
    rc = ''.join([dic[x] for x in sequence])
    return rc

def dna_reverse_complement(sequence):
    return dna_complement(sequence)[::-1]