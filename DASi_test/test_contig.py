'''
Project: jdna
File: test_contig
Author: Justin
Date: 2/23/17

Description: 

'''
import itertools

import numpy as np
import pytest

from DAS.das_contig import Contig, ContigContainer, ContigError
from DASi_test.test_assembly import contig_example


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
            with pytest.raises(ContigError):
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

    with pytest.raises(ContigError):
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

def test_default_contig():