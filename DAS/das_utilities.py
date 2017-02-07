'''
Project: jdna
File: das_utilities
Author: Justin
Date: 2/6/17

Description: 

'''

from das_seqio import *
import numpy as np

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

# def fill_in_gaps(alignments):
#     pairs = itertools.permutations(alignments, 2)
#     gaps = []
#     for l, r in pairs:
#         l_end = l['q_end']
#         r_start = r['q_start']
#         if l_end < r_start:
#             start = l_end + 1
#             end = r_start - 1
#             if start <= end:
#                 g = {'q_start': start, 'q_end': end, 'align_type': 'gap'}
#                 assign_id(g)
#                 gaps.append(g)
#     print '{} gaps'.format(len(gaps))
#     alignments += sorted(gaps, key=lambda x: x['q_start'])

