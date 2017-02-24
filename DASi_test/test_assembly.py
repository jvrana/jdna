'''
Project: jdna
File: test
Author: Justin
Date: 2/16/17

Description: 

'''

from DAS.das_assembly import *
from DAS.das_blast import *
from test_contig import contig_example

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
    contig1 = contig_example.copy()
    contig2 = contig_example.copy()

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