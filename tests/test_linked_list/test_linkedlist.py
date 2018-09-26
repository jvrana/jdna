from copy import copy

import pytest

from jdna.core import DoubleLinkedList


def test_DoubleLinkedList():
    seq = 'This is a test string for the linked data set.'
    l = DoubleLinkedList(sequence=seq)
    assert str(l) == seq
    assert l.head, 'T'
    assert l.head.find_last(), '.'
    assert len(l) == len(seq)


def test_indexing():
    seq = 'This is a test string for the linked data set.'
    l = DoubleLinkedList(sequence=seq)

    assert l[0] == "T"
    assert l[-1] == "."
    assert l[10] == seq[10]


def test_slicing():
    seq = 'This is a test string for the linked data set.'
    l = DoubleLinkedList(sequence=seq)

    assert seq[1:10] == str(l[1:10])


def test_slicing_cyclic():
    pass


def test_cyclic_vs_liner():
    seq = 'This is a test string for the linked data set.'
    l = DoubleLinkedList(sequence=seq)
    assert ~l.is_cyclic()

    l.make_cyclic()
    assert l.is_cyclic()

    assert str(l) == seq
    assert l.head.data, 'T'
    assert l.head.find_last().data, '.'
    assert len(l) == len(seq)

    l.linearize(1)
    assert ~l.is_cyclic()


def test_reindexing():
    seq = 'This is a test string for the linked data set.'
    l = DoubleLinkedList(sequence=seq)
    l.make_cyclic()

    new_index = 10

    l.linearize(new_index)
    assert ~l.is_cyclic()
    assert str(l) == seq[new_index:] + seq[:new_index]


def test_copy():
    seq = 'This is a test string for the linked data set.'
    l = DoubleLinkedList(sequence=seq)
    l_copy = copy(l)
    assert str(l) == str(l_copy)
    assert id(l_copy.links[0]) is not id(l.links[0])
    from copy import deepcopy
    with pytest.raises(NotImplementedError):
        deepcopy(l)


def test_insertion():
    seq = 'This is a test string for the linked data set.'
    l = DoubleLinkedList(sequence=seq)
    for i in range(len(l) + 1):
        l_copy = copy(l)
        insertion_seq = 'XYZ'
        insertion = DoubleLinkedList(sequence=insertion_seq)
        l_copy.insert(insertion, i)
        assert str(l_copy) == seq[:i] + str(insertion_seq) + seq[i:]


def test_insertion_index_error():
    seq = 'This is a test string for the linked data set.'
    l = DoubleLinkedList(sequence=seq)
    with pytest.raises(IndexError):
        l.insert(DoubleLinkedList(sequence='XYZ'), -1)

    with pytest.raises(IndexError):
        l.insert(DoubleLinkedList(sequence='XYZ'), len(l) + 1)


def test_circular_cutting():
    seq = 'This is a test string for the linked data set.'
    l = DoubleLinkedList(sequence=seq)
    l.make_cyclic()
    for cs in [[10, 15, 20], [4, 10, 15]]:
        fragments = [str(x) for x in l.cut(cs)]
        expected = [seq[cs[-1]:] + seq[:cs[0]]]
        expected += [seq[cs[i]:cs[i + 1]] for i in range(len(cs) - 1)]
        expected.sort()
        fragments.sort()
        assert len(expected) == len(fragments)
        for e, f in zip(expected, fragments):
            assert str(e) == str(f)


def test_linear_cutting():
    seq = 'This is a test string for the linked data set.'
    l = DoubleLinkedList(sequence=seq)
    for cs in [[10, 15, 20], [4, 10, 15]]:
        fragments = [str(x) for x in l.cut(cs)]
        expected = [seq[:cs[0]]]
        expected += [seq[cs[i]:cs[i + 1]] for i in range(len(cs) - 1)]
        expected += [seq[cs[-1]:]]
        assert len(expected) == len(fragments)
        expected.sort()
        fragments.sort()
        for e, f in zip(expected, fragments):
            assert str(e) == str(f)
            assert len(fragments) == len(cs) + 1


def test_search():
    query_seq = 'ABCDEFGHIJKLMNOP'
    template = DoubleLinkedList(sequence=query_seq)
    for i in range(len(template) + 1):
        for j in range(i + 1, len(template) + 1):
            query = DoubleLinkedList(sequence=query_seq[i:j])
            assert (i, template.links[i]) == template.search_all(DoubleLinkedList(sequence=query_seq[i:j]))[0]

    query = DoubleLinkedList(sequence='ABCDFG')
    assert [] == template.search_all(query)

    template = DoubleLinkedList(sequence='XXXXGHHHXHGG')
    query = DoubleLinkedList(sequence='XX')
    assert 3 == len(template.search_all(query))

    template = DoubleLinkedList(sequence='longer sentence now for this sentence to find the sentence sentence')
    query = DoubleLinkedList(sequence='sentence')
    assert 4 == len(template.search_all(query))


def test_circular_search():
    # Search Circular
    template = DoubleLinkedList(sequence='XXXXGHHHXHGG')
    template.make_cyclic()
    query = DoubleLinkedList(sequence='HGGXX')
    assert 1 == len(template.search_all(query))


def test_reverse():
    seq = 'XXXXGHHHXHGG'
    l = DoubleLinkedList(sequence=seq)
    l.reverse()
    assert str(l) == seq[::-1]

    seq = 'XXXXGHHHXHGG'
    l = DoubleLinkedList(sequence=seq)
    l.make_cyclic()
    l.reverse()
    assert str(l) == seq[::-1]

    seq = 'XXXXGHHHXHGG'
    l = DoubleLinkedList(sequence=seq)
    l.reverse()
    l.make_cyclic()
    assert str(l) == seq[::-1]


def test_remove():
    seq = 'XXXXGHHHXHGG'
    l = DoubleLinkedList(sequence=seq)
    for i in range(len(seq)):
        l_copy = copy(l)
        l_copy.remove(i)
        if i == 0:
            assert str(l_copy) == seq[i + 1:]
        elif i == len(seq):
            assert str(l_copy) == seq[:-1]
        else:
            assert str(l_copy) == seq[:i] + seq[i + 1:]
