### Test linked set
from nose.tools import *
from hydradna.tests import test_anneal, test_linkedset, test_features, test_reaction, test_convert
from hydradna import *

# print "Testing HydraDNA"
#
# # Test feature creation
# template = Sequence(sequence='AGTGCAGTCTCGGCTATTGTGTCTTGTGATGTCTCTGTTAGTGTCGTA')
# new_feature = Feature('new_feature')
# template.add_feature(5, 10, new_feature)
# print template.get_features()
# frags = template.cut(7)
# for f in frags:
#     print f.get_features()
#
# # Test
# template = Sequence(sequence='AGTGCAGTCTCGGCTATTGTGTCTTGTGATGTCTCTGTTAGTGTCGTA')
# new_feature = Feature('new_feature')
# template.add_feature(5, 10, new_feature)
# template.reverse_complement()
# print template
# print template.get_features()
# frags = template.cut(len(template)-7)
# for f in frags:
#     print f.get_features()


test_anneal.test_anneal()
test_linkedset.test_linkedset()
#test_features.test_features()
test_reaction.test_reaction()
test_convert.test_convert()

# print "Importing benchling api..."
# from benchlingapi import BenchlingAPI
# from benchlingapi.convert import *
#
# print "Uploading Benchling API..."
# api = BenchlingAPI('sk_GbNYfhnukDU30J5fAebIjEj0d4YlJ')
# seq = api.find_sequence('RGR', regex=True)
# print seq
