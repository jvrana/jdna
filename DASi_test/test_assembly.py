'''
Project: jdna
File: test
Author: Justin
Date: 2/16/17

Description: 

'''

from DAS.das_assembly import *
from DAS.das_blast import *

contig_example = {
          "s_end": 757,
          "identical": 155,
          "query_acc": "pINS-011-pEF1a-hcsy4-T",
          "subject_acc": "pRIAS_(CC#15)",
          "query_seq": "CAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCT",
          "contig_type": "blast",
          "start_label": "new_primer",
          "gaps": 0,
          "end_label": "new_primer",
          "bit_score": 310,
          "quality": 1.0,
          "subject_strand": "minus",
            "subject_circular": True,
            "query_circular": True,
          "q_start": 1,
          "q_end": 155,
          "s_start": 911,
          "evalue": 2.95e-84,
          "gap_opens": 0,
          "filename": "templates/pRIAS (CC15).gb",
          "subject_seq": "CAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCT",
          "score": 155,
          "contig_id": 268,
          "subject_length": 9795,
          "query_length": 22240,
          "alignment_length": 155,
          "circular": True
        }

# TODO: fix this reverse stuff...



# TODO: Reverse regions are confusing, fix this


def test_pcr_products():
    pass

def test_expand():
    pass

def test_rc():
    seq =   'agtcGaTcgaN'
    c_seq = 'tcagCtAgctN'
    rc_seq = c_seq[::-1]
    assert c_seq == dna_complement(seq)
    assert rc_seq == dna_reverse_complement(seq)

def test_assembly_graph():
    contig1 = Contig(**contig_example)
    contig2 = Contig(**contig_example)

    contig1.query._Region__length = 10000
    contig2.query._Region__length = contig1.query.length
    contig1.query.start = 100
    contig1.query.end = 200
    contig2.query.start = 201
    contig2.query.end = 300

    print contig1.query.length, contig1.query.start, contig1.query.end

    assert Assembly.assembly_condition(contig1, contig2)

    contig1.query._Region__length = 10000
    contig2.query._Region__length = contig1.query.length
    contig1.query.start = 100
    contig1.query.end = 200
    contig2.query.start = 190
    contig2.query.end = 300

    print contig1.query.length, contig1.query.start, contig1.query.end

    assert Assembly.assembly_condition(contig1, contig2)

    contig1.query._Region__length = 10000
    contig2.query._Region__length = contig1.query.length
    contig1.query.start = 100
    contig1.query.end = 200
    contig2.query.start = 240
    contig2.query.end = 300

    print contig1.query.length, contig1.query.start, contig1.query.end

    assert Assembly.assembly_condition(contig1, contig2)

    contig1.query._Region__length = 1000
    contig2.query._Region__length = contig1.query.length
    contig1.query.start = 900
    contig1.query.end = 950
    contig2.query.start = 10
    contig2.query.end = 100

    print contig1.query.length, contig1.query.start, contig1.query.end

    assert not Assembly.assembly_condition(contig1, contig2)