import random
from copy import copy

import pytest

from jdna import Sequence


@pytest.fixture(scope='function')
def test_str():
    return '^This is a test string for the sequence data set.'


@pytest.fixture(scope='function')
def l(test_str):
    return Sequence(test_str)


# <editor-fold desc="Basic Tests">
def test_Sequence(test_str, l):
    assert str(l) == test_str
    assert l.head.data, 'T'
    assert l.head.find_last().data, '.'
    assert len(l) == len(test_str)


def test_cyclic_vs_liner(test_str, l):
    assert ~l.cyclic

    l.circularize()
    assert l.cyclic

    assert str(l) == test_str
    assert l.head.data, 'T'
    assert l.head.find_last().data, '.'
    assert len(l) == len(test_str)

    l.linearize(1)
    assert ~l.cyclic


def test_reindexing(test_str, l):
    l.circularize()

    new_index = 10

    l.linearize(new_index)
    assert str(l) == test_str[new_index:] + test_str[:new_index]
    assert not l.cyclic


def test_copy(test_str, l):
    l.name = 'sequence name'
    l_copy = copy(l)
    assert str(l) == str(l_copy)
    assert l.name == l_copy.name
    assert l_copy.nodes[0] is not l.nodes[0]
    for n1, n2 in zip(l.nodes, l_copy.nodes):
        assert n1 is not n2
    assert l.__class__ == l_copy.__class__
    from copy import deepcopy
    with pytest.raises(NotImplementedError):
        deepcopy(l)


def test_insertion(test_str, l):
    for i in range(len(l) + 1):
        l_copy = copy(l)
        insertion = Sequence(test_str)
        l_copy.insert(insertion, i)
        assert str(l_copy) == test_str[:i] + str(insertion) + test_str[i:]


def test_insertion2(test_str, l):
    l2 = Sequence(test_str)

    for i in range(len(test_str)):
        assert str(copy(l).insert(copy(l2), i)) == str(test_str)[:i] + str(l2) + str(l)[i:]


def test_insertion_index_error(test_str, l):
    with pytest.raises(IndexError):
        l.insert(Sequence(test_str), -1)

    with pytest.raises(IndexError):
        l.insert(Sequence(test_str), len(l) + 1)


def test_circular_cutting(test_str, l):
    l.circularize()
    for cs in [[10, 15, 20], [4, 10, 15]]:
        fragments = [str(x) for x in l.cut(cs)]
        expected = [test_str[cs[-1]:] + test_str[:cs[0]]]
        expected += [test_str[cs[i]:cs[i + 1]] for i in range(len(cs) - 1)]
        assert set(expected) == set(fragments)


def test_linear_cutting(test_str, l):
    for cs in [[10, 15, 20], [4, 10, 15]]:
        fragments = [str(x) for x in l.cut(cs)]
        expected = [test_str[:cs[0]]]
        expected += [test_str[cs[i]:cs[i + 1]] for i in range(len(cs) - 1)]
        expected += [test_str[cs[-1]:]]
        expected = set(expected)
        fragments = set(fragments)
        assert expected == fragments


def test_reverse():
    seq_str = 'XXXXGHHHXHGG'
    l = Sequence(seq_str)
    l.reverse()
    assert str(l) == seq_str[::-1]

    seq_str = 'XXXXGHHHXHGG'
    l = Sequence(seq_str)
    l.circularize()
    l.reverse()
    assert str(l) == seq_str[::-1]

    seq_str = 'XXXXGHHHXHGG'
    l = Sequence(seq_str)
    l.reverse()
    l.circularize()
    assert str(l) == seq_str[::-1]


# </editor-fold>


"""
feature cutting
copying features
replacing features
copy features from

"""


# <editor-fold desc="DNA Manipulation methods">
def test_cutting():
    seq_str = 'agtcgtatgctgcggcgattctgatgctgatgctgatgtcggta'
    seq = Sequence(seq_str)

    def test_cut(c):
        expected = [seq_str[:c[0]]]
        prev = c[0]
        for loc in c[1:]:
            expected.append(seq_str[prev:loc])
            prev = loc
        expected += [seq_str[c[-1]:]]
        while '' in expected:
            expected.remove('')
        fragments = seq.cut(c)
        assert set([str(x) for x in fragments]) == set(expected)

    test_cut([1, 5, 7])
    test_cut((5, len(seq) - 1))
    test_cut([0, 1, 3, len(seq) - 1, len(seq)])
    test_cut([0, 5, 7])
    test_cut([1, 8, len(seq)])
    assert str(seq) == seq_str


def test_cutting_next():
    seq_str = 'agtcgtatgctgcggcgattctgatgctgatgctgatgtcggta'
    seq = Sequence(seq_str)

    def test_cut(c):
        expected = [seq_str[:c[0] + 1]]
        prev = c[0]
        for loc in c[1:]:
            expected.append(seq_str[prev + 1:loc + 1])
            prev = loc
        expected += [seq_str[c[-1] + 1:]]
        while '' in expected:
            expected.remove('')
        fragments = seq.cut(c, cut_prev=False)
        assert [str(x) for x in fragments] == expected

    test_cut([1, 5, 7])
    test_cut((5, len(seq) - 1))
    with pytest.raises(IndexError):
        test_cut([0, 1, 3, len(seq) - 1, len(seq)])
    test_cut([0, 5, 7])
    test_cut([1, 8, len(seq) - 1])
    assert str(seq) == seq_str


def test_cyclic_cutting():
    seq_str = 'agtcgtatgctgcggcgattctgatgctgatgctgatgtcggta'
    seq = Sequence(seq_str)
    seq.circularize()

    def test_cut(c):
        expected = [seq_str[c[-1]:] + seq_str[:c[0]]]
        prev = c[0]
        for loc in c[1:]:
            expected.append(seq_str[prev:loc])
            prev = loc
        while '' in expected:
            expected.remove('')
        fragments = seq.cut(c)
        for f in fragments:
            assert ~f.cyclic
        assert set([str(x) for x in fragments]) == set(expected)

    test_cut([1])
    test_cut([1, 5, 7])
    test_cut((5, len(seq) - 1))
    test_cut([0, 1, 3, len(seq) - 1, len(seq)])
    test_cut([0, 5, 7])
    test_cut([1, 8, len(seq)])
    assert str(seq) == seq_str


def test_complement():
    seq_str = 'AGTCAGTC'
    seq = Sequence(seq_str)
    seq.complement()
    assert str(seq) == 'TCAGTCAG'
    seq.complement()
    assert str(seq) == seq_str


def test_cyclic_complement():
    seq_str = 'AGTcaGTC'
    seq = Sequence(seq_str)
    seq.circularize()
    seq.complement()
    assert str(seq) == 'TCAgtCAG'


def test_reverse_complement():
    seq_str = 'AGTCAGTC'
    seq = Sequence(seq_str)
    seq.reverse_complement()
    assert str(seq) == 'TCAGTCAG'[::-1]
    seq.reverse_complement()
    assert str(seq) == seq_str


def test_cyclic_reverse_complement():
    seq_str = 'AGTCAGTC'
    seq = Sequence(seq_str)
    seq.circularize()
    seq.reverse_complement()
    assert str(seq) == 'TCAGTCAG'[::-1]


def test_chop(test_str):
    seq = Sequence(test_str)
    for i in range(len(seq)):
        assert str(seq.chop_off_fiveprime(i)) == str(seq)[i:]

    for i in range(len(seq)):
        assert str(seq.chop_off_threeprime(i)) == str(seq)[:i + 1]


@pytest.mark.parametrize('length', [
    0,
    1,
    10,
    100,
    2000
])
def test_random(length):
    if length == 0:
        assert Sequence.random(length) is None
        return
    seq = Sequence.random(length)
    assert len(seq) == length
