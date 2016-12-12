# from nose.tools import *
#
# from hydradna import *
#
#
#
#
# def test_features():
#     # Sequence
#     sequence = 'AGCTGCGGCGGTAGTTTGCTTGAG'
#     s = Sequence(sequence=sequence)
#     s.make_cyclic()
#     assert_true(s.is_cyclic())
#     s.linearize(0)
#     assert_false(s.is_cyclic())
#     s.make_cyclic()
#
#     s_copy = deepcopy(s)
#     s_copy.make_cyclic()
#     s_copy.reverse()
#     s_copy.reverse_complement()
#     s_copy.complement()
#
#
#     # Complement
#     complement = [Nucleotide.base_pairing[x] for x in sequence]
#     s.complement()
#     assert_true(str(s), complement)
#     s.reverse()
#     assert_true(str(s), complement[::-1])
#     s.complement()
#     s.reverse()
#     assert_true(sequence, str(s))
#
#     def capture_features(seq, feature):
#         pos = seq.get_features()[feature][0]
#         assert_equal(len(pos), 1)
#         return pos[0]
#
#
#     # Standard Features
#     assert_equal(len(s.get_features()), 0)
#
#     feature1 = Feature('first_feature', 'test type')
#     start, end = 4, 7
#     s.add_feature(start, end, feature1)
#     assert_equal(capture_features(s, feature1), list((start,end)))
#
#     feature2 = Feature('second_feature', 'test type')
#     start, end = 5, 8
#     s.add_feature(start, end, feature2)
#     assert_equal(capture_features(s, feature2), list((start, end)))
#
#
#     # Find Features
#     assert_equal(s.find_feature('second_feature'), [feature2])
#     assert_equal(s.find_feature('scond_f'), [])
#
#
#     # Circular Feature
#     start = 8
#     end = 5
#     new_feature = s.create_feature('circular_feature', 'test type', start, end)
#     assert_true(list(range(end+1)) + list(range(start, len(s))), s.get_features()[new_feature])
#
#     # Cut feature next
#     s_original = Sequence(sequence="ATGCGTGGCGGATATATCTCTCT")
#
#     f_start = 0
#     f_end = 10
#     for i in range(1, len(s_original)-1):
#         for j in range(i+1, len(s_original)-1):
#             s = deepcopy(s_original)
#             f = Feature('to be cut', 'test')
#             s.add_feature(f_start, f_end, f)
#             cut_start = i
#             cut_end = j
#             fragments = s.cut((cut_start, cut_end), cut_prev=False)
#
#             def get_endpoints(fragment, feature_name):
#                 features = fragment.find_feature('to_be_cut')
#                 ranges = []
#                 for feature in features:
#                     feature = features[0]
#                     start, end = fragment.get_features()[feature][1]
#                     ranges.append((start,end))
#                 return ranges
#             print get_endpoints(fragments[0], 'to_be_cut')
#             print get_endpoints(fragments[1], 'to_be_cut')
#             print get_endpoints(fragments[2], 'to_be_cut')
# #             l1 = [(x.start, x.end) for x in fragments[0].find_feature('to be cut')]
# #             l2 = [(x.start, x.end) for x in fragments[1].find_feature('to be cut')]
# #             l3 = [(x.start, x.end) for x in fragments[2].find_feature('to be cut')]
# #             x1 = [0, cut_start]
# #             x2 = [cut_start + 1, cut_end]
# #             x3 = [cut_end + 1, None]
# #             for X in [x1, x2, x3]:
# #                 if X[1] >= f_end - f_start:
# #                     X[1] = None
# #                 if X[0] > f_end - f_start:
# #                     X[0] = None
# #                 if X == [None, None]:
# #                     X.remove(None)
# #                     X.remove(None)
# #             x1, x2, x3 = [tuple(x) for x in x1, x2, x3]
# #             l1 = [[]] if l1 == [] else l1
# #             l2 = [[]] if l2 == [] else l2
# #             l3 = [[]] if l3 == [] else l3
# #             x1 = [] if x1 == () else x1
# #             x2 = [] if x2 == () else x2
# #             x3 = [] if x3 == () else x3
# #
# #             assert_equal(x1, l1[0])
# #             assert_equal(x2, l2[0])
# #             assert_equal(x3, l3[0])
# #
# # # Cut feature next
# #     s_original = Sequence(sequence="ATGCGTGGCGGATATATCTCTCT")
# #
# #     f_start = 0
# #     f_end = 10
# #     for i in range(1, len(s_original)-1):
# #         for j in range(i+1, len(s_original)-1):
# #             s = deepcopy(s_original)
# #             f = Feature('to be cut', 'test')
# #             s.add_feature(f_start, f_end, f)
# #             cut_start = i
# #             cut_end = j
# #             fragments = s.cut((cut_start, cut_end), cut_prev=True)
# #             l1 = [(x.start, x.end) for x in fragments[0].find_feature('to be cut')]
# #             l2 = [(x.start, x.end) for x in fragments[1].find_feature('to be cut')]
# #             l3 = [(x.start, x.end) for x in fragments[2].find_feature('to be cut')]
# #             x1 = [0, cut_start-1]
# #             x2 = [cut_start, cut_end-1]
# #             x3 = [cut_end, None]
# #             for X in [x1, x2, x3]:
# #                 if X[1] >= f_end - f_start:
# #                     X[1] = None
# #                 if X[0] > f_end - f_start:
# #                     X[0] = None
# #                 if X == [None, None]:
# #                     X.remove(None)
# #                     X.remove(None)
# #             x1, x2, x3 = [tuple(x) for x in x1, x2, x3]
# #             l1 = [[]] if l1 == [] else l1
# #             l2 = [[]] if l2 == [] else l2
# #             l3 = [[]] if l3 == [] else l3
# #             x1 = [] if x1 == () else x1
# #             x2 = [] if x2 == () else x2
# #             x3 = [] if x3 == () else x3
# #
# #             assert_equal(x1, l1[0])
# #             assert_equal(x2, l2[0])
# #             assert_equal(x3, l3[0])