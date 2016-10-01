from jdna.core import Sequence, Feature, Nucleotide, Reaction, Convert, DoubleLinkedList
from copy import copy, deepcopy
import random
import json
import re
# bseq = None
# with open('test_data/seq1.json') as f:
#     bseq = json.load(f)
# hseq = Convert.from_benchling(bseq)
# print len(hseq)
# #
#
#
# p1 = Sequence(sequence=str(hseq)[10:30])
# p2 = Sequence(sequence=str(hseq)[3000:3020]).reverse_complement()
#
# Reaction.pcr(hseq, p1, p2)

import glob

dnas = []
for file in glob.glob('test_data/Fragment*json'):
    with open(file, 'r') as f:
        bseq = json.load(f)
        dnas.append(Convert.from_benchling(bseq))

assert len(Reaction.cyclic_assembly(dnas)) == 1
# print 'cut'
# import itertools
#

# hseq = Sequence(sequence='a'*10000)
# hseq.cut([1,4,5])


# import itertools
# p1 = Sequence(sequence=str(hseq)[0:20])
# p2 = Sequence(sequence=str(hseq)[2000:2020]).reverse_complement()
#
# ann1 = Reaction.anneal_primer(hseq, p1)
# ann2 = Reaction.anneal_primer(hseq, p2)
#
# f = ann1['F'] + ann2['F']
# r = ann1['R'] + ann2['R']
# print r
# pairs = itertools.product(f, r)
# for pair in pairs:
#     print pair
#     copied = copy(hseq)
#     print 'copied'
#     start = pair[0]['pos'] - pair[0]['len']
#     end = len(hseq) - pair[1]['pos'] + pair[1]['len']
#     print 'cutting'
#     print hseq.cut((start, end))



# ignore anything past len(template)
# for each position call longest match


# match = re.search(str(p1), str(hseq))
# print len(str(p1))
# seq_str = str(hseq)
# rotated = seq_str[-len(p1):] + seq_str[:len(p1)]
# print rotated
# match = re.search(str(p1), rotated)
# print ann1
# print ann2
#
# f = ann1['F'] + ann2['F']
# r = ann1['R'] + ann2['R']
#
# print f
# print r


# def test_long_copy():
#     length = 1000000
#     seq_str = ''.join([random.choice('atcg') for x in range(length)])
#     seq = Sequence(sequence=seq_str)
#     print 'copying'
#     copy(seq)
#     # print 'copying'
#     # def copy_sequence():
#     #     seq.create_feature('new feature', 'type', 0, 100)
#     #     seq_copy = Sequence(sequence=str(seq))
#     #     #features = seq.get_features()
#     #     # for f in features:
#     #     #     position, span = copy(features[f])
#     #     #     seq_copy.add_feature(position[0], position[1], copy(f), span[0])
#     #     return seq_copy
#     # copy_sequence()
#     #
#     # def copy_sequence2():
#     #     return copy(seq)
#     # copy_sequence()
#     # # copy_sequence2()
#
#
# test_long_copy()