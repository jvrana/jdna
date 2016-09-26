import pytest

from nose.tools import *
from jdna.core import Link, DoubleLinkedList
from copy import copy


def test_DoubleLinkedList():
    seq = 'This is a test string for the linked data set.'
    l = DoubleLinkedList(data_sequence=seq)
    assert_equal(str(l), seq)
    assert_equal(l.get_first().data, 'T')
    assert_equal(l.get_first().find_last().data, '.')
    assert_equal(len(l), len(seq))

def test_cyclic_vs_liner():
    seq = 'This is a test string for the linked data set.'
    l = DoubleLinkedList(data_sequence=seq)
    assert ~l.is_cyclic()

    l.make_cyclic()
    assert l.is_cyclic()

    l.linearize(1)
    assert ~l.is_cyclic()

def test_reindexing():
    seq = 'This is a test string for the linked data set.'
    l = DoubleLinkedList(data_sequence=seq)
    l.make_cyclic()

    new_index = 10

    l.linearize(new_index)
    assert ~l.is_cyclic()
    assert str(l) == seq[new_index:] + seq[:new_index]

def test_copy():
    seq = 'This is a test string for the linked data set.'
    l = DoubleLinkedList(data_sequence=seq)
    l_copy = copy(l)
    assert str(l) == str(l_copy)
    assert l_copy.get()[0] is not l.get()[0]

def test_insertion():
    seq = 'This is a test string for the linked data set.'
    l = DoubleLinkedList(data_sequence=seq)
    for i in range(len(l)+1):
        l_copy = copy(l)
        insertion = DoubleLinkedList(data_sequence='XYZ')
        l_copy.insert(insertion, i)
        assert str(l_copy) == seq[:i] + str(insertion) + seq[i:]

def test_insertion_index_error():
    seq = 'This is a test string for the linked data set.'
    l = DoubleLinkedList(data_sequence=seq)
    with pytest.raises(IndexError):
        l.insert(DoubleLinkedList(data_sequence='XYZ'), -1)

    with pytest.raises(IndexError):
        l.insert(DoubleLinkedList(data_sequence='XYZ'), len(l)+1)

def test_circular_cutting():
    seq = 'This is a test string for the linked data set.'
    l = DoubleLinkedList(data_sequence=seq)
    l.make_cyclic()
    for cs in [[10, 15, 20], [4, 10, 15]]:
        fragments = [str(x) for x in l.cut(cs)]
        expected = [seq[cs[-1]:] + seq[:cs[0]]]
        expected += [ seq[cs[i]:cs[i+1]] for i in range(len(cs)-1) ]
        for e, f in zip(expected, fragments):
            assert str(e) == str(f)

def test_linear_cutting():
    seq = 'This is a test string for the linked data set.'
    l = DoubleLinkedList(data_sequence=seq)
    for cs in [[10, 15, 20], [4, 10, 15]]:
        fragments = [str(x) for x in l.cut(cs)]
        expected = [seq[:cs[0]]]
        expected += [ seq[cs[i]:cs[i+1]] for i in range(len(cs)-1) ]
        expected += [seq[cs[-1]:]]
        for e, f in zip(expected, fragments):
            assert_true(str(e) == str(f))
            assert_equal(len(fragments), len(cs)+1)

def test_search():
    query_seq = 'ABCDEFGHIJKLMNOP'
    template = DoubleLinkedList(data_sequence=query_seq)
    for i in range(len(template)+1):
        for j in range(i+1, len(template)+1):
            query = DoubleLinkedList(data_sequence=query_seq[i:j])
            assert_equal((i, template.get()[i]), template.search_all(DoubleLinkedList(data_sequence=query_seq[i:j]))[0])

    query = DoubleLinkedList(data_sequence='ABCDFG')
    assert_equal([], template.search_all(query))

    template = DoubleLinkedList(data_sequence='XXXXGHHHXHGG')
    query = DoubleLinkedList(data_sequence='XX')
    assert_equal(3, len(template.search_all(query)))

    template = DoubleLinkedList(data_sequence='longer sentence now for this sentence to find the sentence sentence')
    query = DoubleLinkedList(data_sequence='sentence')
    assert_equal(4, len(template.search_all(query)))

def test_circular_search():
    #Search Circular
    template = DoubleLinkedList(data_sequence='XXXXGHHHXHGG')
    template.make_cyclic()
    query = DoubleLinkedList(data_sequence='HGGXX')
    assert_equal(1, len(template.search_all(query)))

def test_reverse():
    seq = 'XXXXGHHHXHGG'
    l = DoubleLinkedList(data_sequence=seq)
    l.reverse()
    assert str(l) == seq[::-1]

    seq = 'XXXXGHHHXHGG'
    l = DoubleLinkedList(data_sequence=seq)
    l.make_cyclic()
    l.reverse()
    assert str(l) == seq[::-1]

    seq = 'XXXXGHHHXHGG'
    l = DoubleLinkedList(data_sequence=seq)
    l.reverse()
    l.make_cyclic()
    assert str(l) == seq[::-1]