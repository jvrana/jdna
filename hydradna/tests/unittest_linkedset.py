### Test linked set

from hydradna.tests import test_anneal, test_linkedset, test_features, test_reaction, test_convert

print "Testing HydraDNA"
test_anneal.test_anneal()
test_linkedset.test_linkedset()
test_features.test_features()
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
