'''
Project: jdna
File: test
Author: Justin
Date: 2/16/17

Description: 

'''

from DAS.das_contig import *
import pytest
import numpy as np

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



def test_contig_within():

    contig1 = Contig(**contig_example)
    contig2 = Contig(**contig_example)

    contig1.q_start = 1000
    contig1.q_end = 2000
    contig2.q_start = 1500
    contig2.q_end = 3000

    assert False == contig1.is_within(contig2)

    contig2.q_end = contig1.q_end
    assert contig2.is_within(contig1)
    assert not contig1.is_within(contig2)
    assert not contig2.is_within(contig1, inclusive=False)


    contig2.q_end = 3000
    contig2.q_start = contig1.q_end
    assert not contig2.is_within(contig1, inclusive=False)
    assert not contig2.is_within(contig1, inclusive=True)
    assert not contig1.is_within(contig2, inclusive=False)
    assert not contig1.is_within(contig2, inclusive=True)

    contig2.q_start = contig1.q_start
    contig2.q_end = 1900
    assert contig2.is_within(contig1)

def test_is_equivalent():
    contig1 = Contig(**contig_example)
    contig2 = Contig(**contig_example)

    contig1.q_start = 1000
    contig1.q_end = 2000
    contig2.q_start = 1500
    contig2.q_end = 3000
    assert not contig1.equivalent_location(contig2)
    assert not contig2.equivalent_location(contig1)

    contig1.q_start = contig2.q_start
    assert not contig1.equivalent_location(contig2)
    assert not contig2.equivalent_location(contig1)

    contig1.q_end = contig2.q_end
    assert contig1.equivalent_location(contig2)
    assert contig2.equivalent_location(contig1)

def test_pos_within():
    contig1 = Contig(**contig_example)
    contig1.q_start = 1000
    contig1.q_end = 2000

    pos = 1500
    assert contig1.q_pos_within(pos)

    pos = 2001
    assert not contig1.q_pos_within(pos)

    for pos in [contig1.q_start, contig1.q_end]:
        assert not contig1.q_pos_within(pos, inclusive=False)
    for pos in [contig1.q_start, contig1.q_end]:
        assert contig1.q_pos_within(pos, inclusive=True)

def test_break_contig():
    contig1 = Contig(**contig_example)
    contig1.q_start = 1000
    contig1.q_end = 2000

    with pytest.raises(ContigError) as e:
        contig1.break_contig(3000, 4000)

    with pytest.raises(ContigError) as e:
        contig1.break_contig(1000, 2001)

    with pytest.raises(ContigError) as e:
        contig1.break_contig(999, 1500)

    x = np.random.randint(contig1.q_start, contig1.q_end, size=20)
    y = np.random.randint(contig1.q_start, contig1.q_end, size=20)

    for s, e in zip(x, y):
        if s > e:
            with pytest.raises(ContigError) as e:
                contig1.break_contig(s, e)
        else:
            n = contig1.break_contig(s,e)
            assert n.q_start == s
            assert n.q_end == e
            assert n.contig_id > contig1.contig_id
            assert n.parent_id == contig1.contig_id

    with pytest.raises(ContigError) as e:
        contig1.break_contig(contig1.q_start - 1, contig1.q_end)

    with pytest.raises(ContigError) as e:
        contig1.break_contig(contig1.q_start, contig1.q_end + 1)

def create_contigs(list_of_start_and_ends):
    contigs = []
    for x, y in list_of_start_and_ends:
        c = Contig(**contig_example)
        c.q_start = x
        c.q_end = y
        contigs.append(c)
    return contigs

def test_contig_container():
    contigs = [
        [100, 200],
        [200, 500],
        [600, 2000],
        [50, 190],
        [210, 550]
    ]

    # Test no removal
    contigs = create_contigs(contigs)
    c = ContigContainer(contigs=contigs)
    l = len(c.contigs)
    c.remove_redundant_contigs(remove_within=True, remove_equivalent=True, no_removal_if_different_ends=True)
    assert len(c.contigs) == l

    # Test within, equivalent
    new_contig = Contig(**contig_example)
    new_contig.q_start = 50
    new_contig.q_end = 190
    c.contigs.append(new_contig)
    l = len(c.contigs)
    c.remove_redundant_contigs(remove_within=True, remove_equivalent=False, no_removal_if_different_ends=True)
    assert len(c.contigs) == l - 1

    # Test equivalent
    new_contig = Contig(**contig_example)
    new_contig.q_start = 50
    new_contig.q_end = 190
    c.contigs.append(new_contig)
    l = len(c.contigs)
    c.remove_redundant_contigs(remove_within=False, remove_equivalent=True, no_removal_if_different_ends=False)
    assert len(c.contigs) == l - 1

    # Test within
    new_contig = Contig(**contig_example)
    new_contig.q_start = 60
    new_contig.q_end = 180
    c.contigs.append(new_contig)
    l = len(c.contigs)
    c.remove_redundant_contigs(remove_within=False, remove_equivalent=True, no_removal_if_different_ends=True)
    assert len(c.contigs) == l
    c.remove_redundant_contigs(remove_within=True, remove_equivalent=True, no_removal_if_different_ends=True)
    assert len(c.contigs) == l - 1
    assert not new_contig in c.contigs

    # Test within
    new_contig = Contig(**contig_example)
    new_contig.q_start = 60
    new_contig.q_end = 180
    c.contigs.append(new_contig)
    l = len(c.contigs)
    new_contig.start_label = 912
    c.remove_redundant_contigs(remove_within=False, remove_equivalent=True, no_removal_if_different_ends=True)
    assert len(c.contigs) == l
    c.remove_redundant_contigs(remove_within=True, remove_equivalent=True, no_removal_if_different_ends=True)
    assert len(c.contigs) == l
    assert new_contig in c.contigs

    c.remove_redundant_contigs(remove_within=True, remove_equivalent=True, no_removal_if_different_ends=False)
    assert len(c.contigs) == l - 1
    assert not new_contig in c.contigs
