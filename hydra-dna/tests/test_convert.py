from benchlingapi import BenchlingAPI
from benchlingapi.convert import *

from hydradna import *

def test_convert():
    print "Fetching Benchling API"
    api = BenchlingAPI('sk_GbNYfhnukDU30J5fAebIjEj0d4YlJ')

    seq = api.get_sequence('seq_CpsEk4nJ')
    seq['name'] = 'COPY'
    seq['folder'] = 'lib_gCtuy36n'
    api.create_sequence('copy', 'agtcagaggtag', True, seq['folder'], description='alkjdlfjkalksdjfa;jd')

    payload = {
        'name': 'name',
        'description': 'descript',
        'bases': 'agtagagc',
        'circular': True,
        'folder': 'lib_gCtuy36n',
    }
    return api._post('sequences/', payload)
    # print seq
    # #hseq = Sequence(sequence=seq['bases'])
    # hseq = Convert.from_benchling(seq)
    # p1 = Sequence(sequence='agctaggcgtagctgatgcgatcgatctatcgtaGTACctgatgagtccgtgaggacgaaacgagtaagctcgtcG')
    # p2 = Sequence(sequence='tgataacggactagccttattttaacttgctatttctagct')
    #
    # p3 = Sequence(sequence='agctagaaatagcaagttaaaataaggctagtccgttatca')
    # p4 = Sequence(sequence='tacgatagatcgatcgcatcagctacgcctagcttgttgcccagccggcgccagcgagg')
    #
    # new = Reaction.pcr(hseq, p1, p2)[0]
    # bseq = Convert.to_benchling_json(new)
    # bseq['folder'] = 'lib_gCtuy36n'
    # print bseq['description']
    # new2 = Reaction.pcr(hseq, p3, p4)[0]
    # bseq2 = Convert.to_benchling_json(new2)
    # bseq2['folder'] = 'lib_gCtuy36n'
    # print new
    # print new2
    # api.create_sequence(**bseq)
    # api.create_sequence(**bseq2)
    #
    # new3 = Reaction.cyclic_assembly([new, new2])[0]
    # bseq3 = Convert.to_benchling_json((new3))
    # bseq3['folder'] = 'lib_gCtuy36n'
    # api.create_sequence(**bseq3)




