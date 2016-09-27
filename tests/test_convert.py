from benchlingapi import BenchlingAPI
from benchlingapi.convert import *
from jdna.core import Convert
import pytest
import json

# print 'benchilng api...'
# api = BenchlingAPI('sk_GbNYfhnukDU30J5fAebIjEj0d4YlJ')
# print 'done'
#
# folder = api.find_folder('HYDRA', regex=True)['id']
#
#
# bseq = api.get_sequence('seq_O8MkeuZ6')
# with open('seq1.json', 'w') as outfile:
#     json.dump(bseq, outfile)

def test_convert():
    # The last thing you were doing, save a bseq as json load it and test the jsons conver produces.
    json.loads()
    raise Exception