'''
Project: jdna
File: test_region
Author: Justin
Date: 2/23/17

Description: 

'''
import numpy as np
import pytest

from DAS.das_region import Region, RegionError


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
    r = Region(2, 2, 3, True, start_index=1)
    assert r.bounds_end == 3

    s = 1
    e = 2
    l = 10
    start_index = 1
    r = Region(s, e, l, circular=True, start_index=start_index)
    assert r.bounds_end == 10
    assert r.bounds_start == start_index
    indices = range(start_index-3, start_index + l - 1 + 3)
    values = np.arange(start_index, l+start_index)
    print values
    assert len(values) == l
    for x in range(start_index-3, start_index + l - 1):
        assert r.translate_pos(x) == values[x - start_index]
    for x in range(start_index + l - 1, start_index + l - 1 + 3):
        assert r.translate_pos(x) == values[x - start_index - len(values)]

    with pytest.raises(RegionError):
        Region(500, 100, 1000, circular=False)

    start_index = 10
    r = Region(100, 500, 1000, True, start_index=start_index)
    assert r.bounds_start == start_index
    assert r.bounds_end == start_index + r.length - 1
    for x in range(-5, 300):
        r = Region(100, x, 200, circular=True, start_index=30)
        assert r.end <= r.bounds_end
    for x in range(-5, 300):
        r = Region(x, 100, 200, circular=True, start_index=30)
        assert r.start >= r.bounds_start


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
    assert not r.within_region(0, inclusive=True)
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