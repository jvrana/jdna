from jdna.core import *
#
# seq = 'XXXXGHHHXHGG'
# l = DoubleLinkedList(data_sequence=seq)
# for i in range(len(seq)):
#     l_copy = copy(l)
#     print 'removing', i
#     print l_copy
#     l_copy.remove(i)
#     print l_copy
#     if i == 0:
#         assert str(l_copy) == seq[i + 1:]
#     elif i == len(seq):
#         assert str(l_copy) == seq[:-1]
#     else:
#         print seq[:i - 1] + seq[i + 1:]
#         assert str(l_copy) == seq[:i] + seq[i + 1:]


# seq = Sequence(sequence='agtcga')
# print seq.__dict__
# print type(seq)
# print type(copy(seq))