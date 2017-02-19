'''
Project: jdna
File: test
Author: Justin
Date: 2/16/17

Description: 

'''

from DAS.das_contig import *
from DAS.das_utilities import *
from DAS.das_blast import *
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



def test_region():
    s = 1
    e = 2
    l = 10
    start_index = 0
    r = Region(s, e, l, circular=True, start_index=start_index)
    indices = range(start_index-3, start_index + l - 1 + 3)
    values = np.arange(start_index, l+start_index)
    print values
    assert len(values) == l
    for x in range(start_index-3, start_index + l - 1):
        assert r.translate_pos(x) == values[x - start_index]
        for x in range(start_index + l - 1, start_index + l - 1 + 3):
            assert r.translate_pos(x) == values[x - len(values)]




    # assert r.translate_pos(-1) == -1 + l
    # assert r.translate_pos(l + start_index - 1 + 2) == start_index + 2 - 1



    # for start_index in range(-2, 2):
    #     indices = np.array(range(start_index - 10, l + 10))
    #     new_indices = indices + start_index
    #     for i, pos in enumerate(indices):
    #         start_index = 0
    #         r = Region(s, e, l, circular=True, start_index=start_index)
    #         x = r.translate_pos(pos)
    #         assert x == new_indices[i]


def test_circular_pos():
    pos = 4
    length = 1000
    start_index = 0
    new = convert_circular_position(pos, length, start_index)
    assert new == pos

    pos = 4
    length = 1000
    start_index = 0
    new = convert_circular_position(length + start_index, length, start_index)
    assert new == start_index

    pos = 4
    length = 1000
    start_index = 1
    new = convert_circular_position(pos, length, start_index)
    assert new == pos

    pos = 4
    length = 1000
    start_index = 1
    new = convert_circular_position(length + start_index, length, start_index)
    assert new == start_index

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

    contig1.break_contig(contig1.q_start, contig1.q_start+1)
    contig1.break_contig(contig1.q_end-1, contig1.q_end)
    x = np.random.randint(contig1.q_start, contig1.q_end, size=20)
    y = np.random.randint(contig1.q_start, contig1.q_end, size=20)

    for s, e in zip(x, y):
        if s > e:
            with pytest.raises(ContigError):
                contig1.break_contig(s, e)
        else:
            n = contig1.break_contig(s,e)
            assert n.q_start == s
            assert n.q_end == e
            assert n.contig_id > contig1.contig_id
            assert n.parent_id == contig1.contig_id
            assert n.subject_length == contig1.subject_length
            assert 0 <= n.s_start <= n.subject_length
            assert 0 <= n.s_end <= n.subject_length

    with pytest.raises(ContigError) as e:
        contig1.break_contig(contig1.q_start - 1, contig1.q_end)

    with pytest.raises(ContigError) as e:
        contig1.break_contig(contig1.q_start, contig1.q_end + 1)

    n = contig1.break_contig(1010, 1010)
    assert n.q_start == 1010
    assert n.q_end == 1010

    with pytest.raises(ContigError):
        contig1.break_contig(contig1.q_start+1, contig1.q_start-1)

    with pytest.raises(ContigError):
        contig1.break_contig(contig1.q_start+10, contig1.q_start+9)

    with pytest.raises(ContigError):
        contig1.break_contig(contig1.q_start-10, contig1.q_start-1)

    with pytest.raises(ContigError):
        contig1.break_contig(contig1.q_start-10, contig1.q_start+1)

    n = contig1.break_contig(contig1.q_start + 100, contig1.q_end - 100, start_label="new", end_label="new2")
    assert n.start_label == "new"
    assert n.end_label == "new2"

    n = contig1.break_contig(contig1.q_start + 100, contig1.q_end - 100, start_label="new", end_label=None)
    assert n.start_label == "new"
    assert n.end_label == Contig.DEFAULT_END

    n = contig1.break_contig(contig1.q_start + 100, contig1.q_end - 100, start_label=None, end_label="new2")
    assert n.start_label == Contig.DEFAULT_END
    assert n.end_label == "new2"

def test_subject_break_contig():
    contig1 = Contig(**contig_example)
    contig1.q_start = 1000
    contig1.q_end = 2000
    contig1.s_start = 100
    contig1.s_end = 2100
    contig1.subject_length = 2000

def create_contigs(list_of_start_and_ends):
    contigs = []
    for x, y in list_of_start_and_ends:
        c = Contig(**contig_example)
        c.q_start = x
        c.q_end = y
        contigs.append(c)
    return contigs

def test_contig_overlaps():
    contig1 = Contig(**contig_example)
    contig2 = Contig(**contig_example)

    contig1.q_start = 1000
    contig1.q_end = 2000
    contig2.q_start = 2001
    contig2.q_end = 3000

    assert not contig1.overlaps(contig2)
    assert not contig2.overlaps(contig1)

    contig2.q_start = contig1.q_end
    assert not contig1.overlaps(contig2, inclusive=False)
    assert not contig2.overlaps(contig1, inclusive=False)
    assert contig1.overlaps(contig2, inclusive=True)
    assert contig2.overlaps(contig1, inclusive=True)

    contig2.q_start = 1999
    assert contig1.overlaps(contig2, inclusive=True)
    assert contig2.overlaps(contig1, inclusive=True)
    assert contig1.overlaps(contig2, inclusive=False)
    assert contig2.overlaps(contig1, inclusive=False)

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


def test_divide_contig():
    contig1 = Contig(**contig_example)

    contig1.q_start = 1000
    contig1.q_end = 5000

    start_points = [
        (1001, 514),
    ]
    end_points = [
        (1200, 516)
    ]

    new_contigs = contig1.divide_contig(startpoints=start_points, endpoints=end_points, )
    assert len(new_contigs) == len(start_points)

    end_points.append((5001, 5))
    new_contigs = contig1.divide_contig(startpoints=start_points, endpoints=end_points, )
    assert len(new_contigs) == 1

    end_points.append((contig1.q_start-1, 5))
    new_contigs = contig1.divide_contig(startpoints=start_points, endpoints=end_points, )
    assert len(new_contigs) == 1

    end_points.append((contig1.q_end-100, 500))
    new_contigs = contig1.divide_contig(startpoints=start_points, endpoints=end_points, )
    assert len(new_contigs) == 2

    start_points.append((contig1.q_start+100, 501))
    new_contigs = contig1.divide_contig(startpoints=start_points, endpoints=end_points, )
    assert len(new_contigs) == 4

    # 3x starts + 2x end
    start_points.append((contig1.q_start + 101, 501))
    new_contigs = contig1.divide_contig(startpoints=start_points, endpoints=end_points, )
    assert len(new_contigs) == 6

    start_points = [
        (1001, 514),
        (contig1.q_start, 100)
    ]
    end_points = [
        (1200, 516),
        (contig1.q_end, None)
    ]
    new_contigs = contig1.divide_contig(start_points, end_points, )
    assert len(new_contigs) == 4
    new_contigs = contig1.divide_contig(start_points, end_points, include_contig=False)
    assert len(new_contigs) == 3

    start_points = [
        (1001, 514),
        (1002, None)
    ]
    end_points = [
        (1200, 516),
        (1300, None)
    ]
    new_contigs = contig1.divide_contig(start_points, end_points, include_contig=False)
    pairs = itertools.product(start_points, end_points)
    for contig, pair in zip(new_contigs, pairs):
        start, start_label = pair[0]
        end, end_label = pair[1]
        if start_label == None:
            start_label = Contig.DEFAULT_END
        if end_label == None:
            end_label = Contig.DEFAULT_END
        assert contig.q_start == start
        assert contig.q_end == end
        assert contig.start_label == start_label
        assert contig.end_label == end_label


def test_break_long_contig():
    pass

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

def test_reindexed_templates():
    b = BLAST('db', 'data/blast_test/reindexed_template', 'data/blast_test/reindexed_design/test.gb', 'data/blast_results', 'data/blast_results/results.out', evalue=10.0, ungapped='',
              penalty=-100, perc_identity=100)
    b.makedbfromdir()
    b.runblast()
    contig_container = b.parse_results(contig_type=Contig.TYPE_BLAST)
    assert len(contig_container.contigs) > 1

    contig_container.fuse_circular_fragments()
    contig_container.remove_redundant_contigs(remove_equivalent=True, remove_within=True,
                                              no_removal_if_different_ends=True)
    assert len(contig_container.contigs) == 1
    c = contig_container.contigs[0]
    assert c.q_start == 1
    assert c.q_end == c.query_length
    assert 1 <= c.s_start <= c.subject_length
    assert 1 <= c.s_end <= c.subject_length



def test_blast_same():
    b = BLAST('db', 'data/blast_test/reindexed_template', 'data/blast_test/reindexed_design/test_reindexed.gb', 'data/blast_results', 'data/blast_results/results.out', evalue=10.0, ungapped='',
              penalty=-100, perc_identity=100)
    b.makedbfromdir()
    b.runblast()
    contig_container = b.parse_results(contig_type=Contig.TYPE_BLAST)
    contig_container.fuse_circular_fragments()
    contig_container.remove_redundant_contigs(remove_within=True)

    contigs = contig_container.contigs
    for c in contigs:
        print c.subject_acc, c.q_start, c.q_end, c.query_length, c.s_start, c.s_end, c.subject_length

    assert len(contig_container.contigs) == 2 # because its circular


def test_pseudo_blast():
    design_path = 'data/blast_test/designs/pmodkan-ho-pact1-z4-er-vpr.gb'
    p = BLAST('primerdb', 'data/blast_test/primers', design_path, '', 'database/primerresults.out')
    p.perfect_matches()