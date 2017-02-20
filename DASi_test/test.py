'''
Project: jdna
File: test
Author: Justin
Date: 2/16/17

Description: 

'''

from DAS.das_contig import *
from DAS.das_assembly import *
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
def test_reverse_region():
    s = 1
    e = 5
    l = 10
    start_index = 0

    r = Region(s, e, l, direction=Region.FORWARD, circular=True, start_index=start_index)

    r = Region(e, s, l, direction=Region.REVERSE, circular=False, start_index=start_index)
    r = Region(e, s, l, direction=Region.REVERSE, circular=True, start_index=start_index)
    print r.start, r.end
    assert r.start == s
    assert r.end == e


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

    with pytest.raises(RegionError):
        Region(500, 100, 1000, circular=False)

    start_index = 10
    r = Region(100, 500, 1000, True, start_index=start_index)
    assert r.bounds_start == start_index
    assert r.bounds_end == start_index + r.length - 1


def test_region_within():
    r = Region(1000, 2000, 4000, circular=True, start_index=1)
    assert r.within_region(1000, inclusive=True)
    assert not r.within_region(1000, inclusive=False)
    assert r.within_region(1500, inclusive=True)
    assert r.within_region(2000, inclusive=True)
    assert not r.within_region(2000, inclusive=False)
    assert not r.within_region(2001, inclusive=True)
    assert not r.within_region(999, inclusive=True)

    r = Region(900, 100, 1000, circular=True, start_index=1)
    assert r.within_region(1, inclusive=True)
    assert r.within_region(10, inclusive=True)
    assert r.within_region(901, inclusive=True)
    assert r.within_region(0, inclusive=True)
    assert r.within_region(900, inclusive=True)
    assert not r.within_region(900, inclusive=False)
    assert not r.within_region(100, inclusive=False)
    assert r.within_region(100, inclusive=True)



def test_set_region():
    r = Region(100, 900, 1000, circular=True, start_index=10)

    start = 100
    span = 50
    r.set_region_by_start_and_span(start, span)
    assert r.region_span == span
    assert r.start == start
    assert r.end == r.start + r.region_span - 1

    with pytest.raises(RegionError):
        r.set_region_by_start_and_span(100, 0)

    r.set_region_by_start_and_span(900, 200)
    print r.start, r.end, r.region_span

# TODO: Reverse regions are confusing, fix this
def test_region_copy():
    r = Region(50, 20, 1000, False, start_index=2, name='name', direction=Region.REVERSE)
    r2 = r.copy()
    assert r.same_context(r2)

def test_region_span():
    # TODO: fix span calculation for circular spans
    def f(span, start, length, start_index):
        r = Region.create(length, True, start_index=start_index).set_region_by_start_and_span(start, span)
        print r.start, r.end, r.bounds_start, r.bounds_end
        print r.region_span
        assert r.region_span == span
        return r

    f(200, 200, 1000, 0)
    f(5, 98, 100, 0)
    f(5, 98, 100, 10)
    f(5, 104, 100, 5)

def test_region_gap():
    r = Region(20, 50, 100, False, start_index=0)
    r2 = Region(60, 70, 100, False, start_index=0)
    g = r.get_gap(r2)
    assert g.start == 51
    assert g.end == 59

    r = Region(20, 50, 100, False, start_index=0)
    r2 = Region(51, 70, 100, False, start_index=0)
    assert r.get_gap(r2) is None

    r = Region(90, 99, 100, False, start_index=1)
    r2 = Region(2, 70, 100, False, start_index=1)
    g = r.get_gap(r2)
    assert g.start == 100
    assert g.end == 1

def test_region_overlap():
    r = Region(10, 50, 100, False, start_index=0)
    r2 = Region(20, 90, 100, False, start_index=0)
    assert r.end_overlaps_with(r2)
    assert not r2.end_overlaps_with(r)

    r = Region(10, 50, 100, False, start_index=0)
    r2 = Region(20, 90, 100, True, start_index=0)
    assert not r.end_overlaps_with(r2)
    assert not r2.end_overlaps_with(r)

    r = Region(10, 50, 100, False, start_index=0)
    r2 = Region(20, 90, 100, False, start_index=1)
    assert not r.end_overlaps_with(r2)
    assert not r2.end_overlaps_with(r)

    r = Region(90, 10, 100, True, start_index=2)
    r2 = Region(5, 20, 100, True, start_index=2)
    assert r.end_overlaps_with(r2)
    assert not r2.end_overlaps_with(r)

    r = Region(50, 10, 100, False, start_index=0, direction=Region.REVERSE)
    r2 = Region(90, 20, 100, False, start_index=0, direction=Region.REVERSE)
    assert r.end_overlaps_with(r2)
    assert not r2.end_overlaps_with(r)

    r = Region(10, 50, 100, False, start_index=0)
    r2 = Region(20, 90, 100, False, start_index=0)
    overlap = r.get_overlap(r2)
    print overlap.start, overlap.end
    assert overlap.same_context(r)
    assert overlap.same_context(r2)
    assert overlap.start == 20
    assert overlap.end == 50
    assert overlap.direction == r.direction
    assert r2.get_overlap(r) is None

    r = Region(90, 10, 100, True, start_index=2)
    r2 = Region(95, 20, 100, True, start_index=2)
    overlap = r.get_overlap(r2)
    assert overlap.start == 95
    assert overlap.end == 10

def test_subregion():
    r = Region(20, 50, 100, False, start_index=0)
    r.sub_region(20, 30)
    s = r.sub_region(30,50)
    assert s.start == 30
    assert s.end == 50
    with pytest.raises(RegionError):
        r.sub_region(19, 30)
    with pytest.raises(RegionError):
        r.sub_region(30, 51)

    r = Region(90, 10, 100, True, start_index=0)
    s = r.sub_region(95, 5)
    s_compare = Region(95, 5, 100, True, start_index=0)
    assert s.start == 95
    assert s.end == 5
    assert s_compare.region_span == s.region_span

def test_region_fuse():

    r = Region(20, 50, 100, False, start_index=0)
    r2 = Region(51, 70, 100, False, start_index=0)
    assert r.consecutive_with(r2)
    assert not r2.consecutive_with(r)

    r = Region(20, 50, 100, False, start_index=0)
    r2 = Region(51, 70, 100, False, start_index=1)
    assert not r.consecutive_with(r2)
    assert not r2.consecutive_with(r)

    r = Region(95, 99, 100, False, start_index=0)
    r2 = Region(0, 70, 100, False, start_index=0)
    assert not r.consecutive_with(r2)
    assert not r2.consecutive_with(r)

    r = Region(95, 99, 100, True, start_index=0)
    r2 = Region(0, 70, 100, True, start_index=0)
    assert r.consecutive_with(r2)
    assert not r2.consecutive_with(r)

    r = Region(95, 100, 100, True, start_index=1)
    r2 = Region(1, 70, 100, True, start_index=1)
    assert r.consecutive_with(r2)
    assert not r2.consecutive_with(r)

    r.fuse(r2)
    assert r.start == 95
    assert r.end == 70
    assert r2.start == 1
    assert r2.end == 70

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

def test_contig_within():

    contig1 = Contig(**contig_example)
    contig2 = Contig(**contig_example)

    contig1.query.start = 1000
    contig1.query.end = 2000
    contig2.query.start = 1500
    contig2.query.end = 3000

    assert False == contig1.is_within(contig2)

    contig2.query.end = contig1.query.end
    assert contig2.is_within(contig1)
    assert not contig1.is_within(contig2)
    assert not contig2.is_within(contig1, inclusive=False)


    contig2.query.end = 3000
    contig2.query.start = contig1.query.end
    assert not contig2.is_within(contig1, inclusive=False)
    assert not contig2.is_within(contig1, inclusive=True)
    assert not contig1.is_within(contig2, inclusive=False)
    assert not contig1.is_within(contig2, inclusive=True)

    contig2.query.start = contig1.query.start
    contig2.query.end = 1900
    assert contig2.is_within(contig1)

def test_is_equivalent():
    contig1 = Contig(**contig_example)
    contig2 = Contig(**contig_example)

    contig1.query.start = 1000
    contig1.query.end = 2000
    contig2.query.start = 1500
    contig2.query.end = 3000
    assert not contig1.equivalent_location(contig2)
    assert not contig2.equivalent_location(contig1)

    contig1.query.start = contig2.query.start
    assert not contig1.equivalent_location(contig2)
    assert not contig2.equivalent_location(contig1)

    contig1.query.end = contig2.query.end
    assert contig1.equivalent_location(contig2)
    assert contig2.equivalent_location(contig1)

def test_pos_within():
    contig1 = Contig(**contig_example)
    contig1.query.start = 1000
    contig1.query.end = 2000

    pos = 1500
    assert contig1.query.within_region(pos)

    pos = 2001
    assert not contig1.query.within_region(pos)

    for pos in [contig1.query.start, contig1.query.end]:
        assert not contig1.query.within_region(pos, inclusive=False)
    for pos in [contig1.query.start, contig1.query.end]:
        assert contig1.query.within_region(pos, inclusive=True)

def test_break_contig():
    contig1 = Contig(**contig_example)
    contig1.query.start = 1000
    contig1.query.end = 2000
    contig1.query._Region__circular = False

    with pytest.raises(ContigError) as e:
        contig1.break_contig(3000, 4000)

    with pytest.raises(ContigError) as e:
        contig1.break_contig(1000, 2001)

    with pytest.raises(ContigError) as e:
        contig1.break_contig(999, 1500)

    contig1.break_contig(contig1.query.start, contig1.query.start+1)
    contig1.break_contig(contig1.query.end-1, contig1.query.end)
    x = np.random.randint(contig1.query.start, contig1.query.end, size=20)
    y = np.random.randint(contig1.query.start, contig1.query.end, size=20)

    for s, e in zip(x, y):
        if s > e:
            with pytest.raises(RegionError):
                contig1.break_contig(s, e)
        else:
            n = contig1.break_contig(s,e)
            assert n.query.start == s
            assert n.query.end == e
            assert n.contig_id > contig1.contig_id
            assert n.parent_id == contig1.contig_id
            assert n.query.length == contig1.query.length
            assert 0 <= n.subject.start <= n.query.length
            assert 0 <= n.subject.end <= n.query.length

    with pytest.raises(ContigError) as e:
        contig1.break_contig(contig1.query.start - 1, contig1.query.end)

    with pytest.raises(ContigError) as e:
        contig1.break_contig(contig1.query.start, contig1.query.end + 1)

    n = contig1.break_contig(1010, 1010)
    assert n.query.start == 1010
    assert n.query.end == 1010

    with pytest.raises(ContigError):
        contig1.break_contig(contig1.query.start+1, contig1.query.start-1)

    with pytest.raises(RegionError):
        contig1.break_contig(contig1.query.start+10, contig1.query.start+9)

    with pytest.raises(ContigError):
        contig1.break_contig(contig1.query.start-10, contig1.query.start-1)

    with pytest.raises(ContigError):
        contig1.break_contig(contig1.query.start-10, contig1.query.start+1)

    n = contig1.break_contig(contig1.query.start + 100, contig1.query.end - 100, start_label="new", end_label="new2")
    assert n.start_label == "new"
    assert n.end_label == "new2"

    n = contig1.break_contig(contig1.query.start + 100, contig1.query.end - 100, start_label="new", end_label=None)
    assert n.start_label == "new"
    assert n.end_label == Contig.DEFAULT_END

    n = contig1.break_contig(contig1.query.start + 100, contig1.query.end - 100, start_label=None, end_label="new2")
    assert n.start_label == Contig.DEFAULT_END
    assert n.end_label == "new2"

def test_subject_break_contig():
    contig1 = Contig(**contig_example)
    contig1.query.start = 1000
    contig1.query.end = 2000
    contig1.subject.start = 100
    contig1.subject.end = 2100
    contig1.query.__length = 4000

def create_contigs(list_of_start_and_ends):
    contigs = []
    for x, y in list_of_start_and_ends:
        c = Contig(**contig_example)
        c.query.start = x
        c.query.end = y
        contigs.append(c)
    return contigs

def test_contig_overlaps():
    contig1 = Contig(**contig_example)
    contig2 = Contig(**contig_example)

    contig1.query.start = 1000
    contig1.query.end = 2000
    contig2.query.start = 2001
    contig2.query.end = 3000

    assert not contig1.overlaps(contig2)
    assert not contig2.overlaps(contig1)

    contig2.query.start = contig1.query.end
    assert not contig1.overlaps(contig2, inclusive=False)
    assert not contig2.overlaps(contig1, inclusive=False)
    assert contig1.overlaps(contig2, inclusive=True)
    assert contig2.overlaps(contig1, inclusive=True)

    contig2.query.start = 1999
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
    new_contig.query.start = 50
    new_contig.query.end = 190
    c.contigs.append(new_contig)
    l = len(c.contigs)
    c.remove_redundant_contigs(remove_within=True, remove_equivalent=False, no_removal_if_different_ends=True)
    assert len(c.contigs) == l - 1

    # Test equivalent
    new_contig = Contig(**contig_example)
    new_contig.query.start = 50
    new_contig.query.end = 190
    c.contigs.append(new_contig)
    l = len(c.contigs)
    c.remove_redundant_contigs(remove_within=False, remove_equivalent=True, no_removal_if_different_ends=False)
    assert len(c.contigs) == l - 1

    # Test within
    new_contig = Contig(**contig_example)
    new_contig.query.start = 60
    new_contig.query.end = 180
    c.contigs.append(new_contig)
    l = len(c.contigs)
    c.remove_redundant_contigs(remove_within=False, remove_equivalent=True, no_removal_if_different_ends=True)
    assert len(c.contigs) == l
    c.remove_redundant_contigs(remove_within=True, remove_equivalent=True, no_removal_if_different_ends=True)
    assert len(c.contigs) == l - 1
    assert not new_contig in c.contigs

    # Test within
    new_contig = Contig(**contig_example)
    new_contig.query.start = 60
    new_contig.query.end = 180
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

    contig1.query.start = 1000
    contig1.query.end = 5000

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

    end_points.append((contig1.query.start-1, 5))
    new_contigs = contig1.divide_contig(startpoints=start_points, endpoints=end_points, )
    assert len(new_contigs) == 1

    end_points.append((contig1.query.end-100, 500))
    new_contigs = contig1.divide_contig(startpoints=start_points, endpoints=end_points, )
    assert len(new_contigs) == 2

    start_points.append((contig1.query.start+100, 501))
    new_contigs = contig1.divide_contig(startpoints=start_points, endpoints=end_points, )
    assert len(new_contigs) == 4

    # 3x starts + 2x end
    start_points.append((contig1.query.start + 101, 501))
    new_contigs = contig1.divide_contig(startpoints=start_points, endpoints=end_points, )
    assert len(new_contigs) == 6

    start_points = [
        (1001, 514),
        (contig1.query.start, 100)
    ]
    end_points = [
        (1200, 516),
        (contig1.query.end, None)
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
        assert contig.query.start == start
        assert contig.query.end == end
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

    contigs = contig_container.contigs
    contig_container.fuse_circular_fragments()
    contig_container.remove_redundant_contigs(remove_equivalent=True, remove_within=True,
                                              no_removal_if_different_ends=True)
    for c in contigs:
        print c.query.name, c.query.start, c.query.end, c.query.length, c.subject.start, c.subject.end, c.query.length
    assert len(contig_container.contigs) == 1
    c = contig_container.contigs[0]
    assert c.subject.region_span == c.query.region_span
    assert c.subject.region_span == c.query.length



def test_blast_same():
    b = BLAST('db', 'data/blast_test/reindexed_template', 'data/blast_test/reindexed_design/test_reindexed.gb', 'data/blast_results', 'data/blast_results/results.out', evalue=10.0, ungapped='',
              penalty=-100, perc_identity=100)
    b.makedbfromdir()
    b.runblast()
    contig_container = b.parse_results(contig_type=Contig.TYPE_BLAST)

    contig_container.fuse_circular_fragments()
    contig_container.remove_redundant_contigs(remove_within=True)
    # contig_container.remove_redundant_contigs(remove_within=True)

    contigs = contig_container.contigs
    for c in contigs:
        print c.query.name, c.query.start, c.query.end, c.query.length, c.subject.start, c.subject.end, c.query.length

    # Fuse the fragments and there should only be one contig
    assert len(contig_container.contigs) == 1

def test_blast_reverse_complement():
    b = BLAST('db', 'data/blast_test/rc_template', 'data/blast_test/rc_design/pnl1105_pgp8g2-ccdb.gb', 'data/blast_results', 'data/blast_results/results.out', evalue=10.0, ungapped='',
              penalty=-100, perc_identity=100)
    b.makedbfromdir()
    b.runblast()
    contig_container = b.parse_results(contig_type=Contig.TYPE_BLAST)

    contig_container.fuse_circular_fragments()
    contig_container.remove_redundant_contigs(remove_within=True)

    contigs = contig_container.contigs
    for c in contigs:
        print c.query.name, c.query.start, c.query.end, c.query.length, c.subject.start, c.subject.end, c.query.length, c.subject.direction



def test_pseudo_blast():
    design_path = 'data/blast_test/designs/pmodkan-ho-pact1-z4-er-vpr.gb'
    p = BLAST('primerdb', 'data/blast_test/primers', design_path, '', 'database/primerresults.out')
    print p.perfect_matches()

def test_primer_binding():
    assert 1 == 2

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