from copy import copy

import pytest

from jdna.linked_list import DoubleLinkedList


def test_empty():
    l = DoubleLinkedList()
    l2 = DoubleLinkedList()
    l3 = l + l2
    for n in l3:
        print(n)


def test_linked_list_constructor(test_str, linked_list):
    assert str(linked_list) == test_str
    assert linked_list.head, test_str[0]
    assert linked_list.head.find_last(), test_str[-1]
    assert len(linked_list) == len(test_str)


def test_indexing(test_str, linked_list):
    assert linked_list[0].data == test_str[0]
    assert linked_list[-1].data == test_str[-1]
    assert linked_list[10].data == test_str[10]


@pytest.mark.parametrize('i,j', [
    (0, 10),
    (2, 5),
    (5, 2),
    (1, 1),
    (1, 2)
])
def test_slicing(i, j, test_str, linked_list):
    mylist = linked_list[i:j]
    if i >= j:
        assert mylist.is_empty()
    else:
        assert str(mylist) == test_str[i:j]


@pytest.mark.parametrize('i,j', [
    (-5, 5),
    (5, 2),
    (1, 1),
])
def test_slicing_cyclic(i, j, test_str, linked_list):
    linked_list.cyclic = True
    if i == j:
        assert linked_list[i:j].is_empty()
    else:
        assert str(linked_list[i:j]) == test_str[i:] + test_str[:j]


def test_cyclic_vs_linear(test_str, linked_list):
    assert ~linked_list.cyclic

    linked_list.circularize()
    assert linked_list.cyclic

    assert str(linked_list) == test_str
    assert linked_list.head.data, 'T'
    assert linked_list.head.find_last().data, '.'
    assert len(linked_list) == len(test_str)

    linked_list.linearize(1)
    assert ~linked_list.cyclic


@pytest.mark.parametrize('index', [
    0,
    5,
    10
])
def test_reindexing(index, test_str, linked_list):
    linked_list.circularize()
    linked_list.linearize(index)
    assert ~linked_list.cyclic
    assert str(linked_list) == test_str[index:] + test_str[:index]


def test_copy(test_str, linked_list):
    l_copy = copy(linked_list)
    assert str(linked_list) == str(l_copy)
    assert id(l_copy.nodes[0]) is not id(linked_list.nodes[0])
    from copy import deepcopy
    with pytest.raises(NotImplementedError):
        deepcopy(linked_list)


@pytest.mark.parametrize('i', [None, 1, 4, 6, 8])
@pytest.mark.parametrize('j', [None, 1, 4, 6, 8])
def test_new_slice(test_str, linked_list, i, j):
    n1 = linked_list.get(i)
    n2 = linked_list.get(j)

    new_slice = DoubleLinkedList.new_slice(n1, n2)
    if i is None and j is None:
        assert new_slice.is_empty()
    elif i is None:
        assert str(new_slice) == str(linked_list)[:j+1]
    elif j is None:
        assert str(new_slice) == str(linked_list)[i:]
    elif i > j:
        assert new_slice.is_empty()
    else:
        assert str(new_slice) == str(linked_list)[i:j+1]


def test_new_slice_with_node_never_reached(test_str, linked_list):
    n1 = linked_list.get(1)
    n2 = linked_list.copy().get(3)
    new_slice = DoubleLinkedList.new_slice(n1, n2)
    assert new_slice.is_empty()

def test_insertion(test_str, linked_list):
    for i in range(len(linked_list) + 1):
        l_copy = copy(linked_list)
        insertion_test_str = 'XYZ'
        insertion = DoubleLinkedList(data=insertion_test_str)
        l_copy.insert(insertion, i)
        assert str(l_copy) == test_str[:i] + str(insertion_test_str) + test_str[i:]


def test_insertion_index_error(test_str, linked_list):
    with pytest.raises(IndexError):
        linked_list.insert(DoubleLinkedList(data='XYZ'), -1)

    with pytest.raises(IndexError):
        linked_list.insert(DoubleLinkedList(data='XYZ'), len(linked_list) + 1)


@pytest.mark.parametrize('cut_sites,circular,expected_num_fragments', [
    ([10], False, 2),
    ([10], True, 1),
    ([10, 15, 20], False, 4),
    ([10, 15, 20], True, 3),
    ([4, 10, 15], False, 4),
    ([4, 10, 15], True, 3),
])
def test_cutting(cut_sites, circular, expected_num_fragments, test_str):
    linkedlist = DoubleLinkedList(data=test_str)
    linkedlist.cyclic = circular
    assert linkedlist.cyclic == circular, "fragment should be cyclic={}".format(circular)
    fragments = linkedlist.cut(cut_sites)
    assert len(fragments) == expected_num_fragments, "Wrong number of cut fragments"
    for f in fragments:
        assert not f.cyclic, "all cut fragments should be linear"

    expected = []
    expected.append(test_str[:cut_sites[0]])
    if len(cut_sites) > 1:
        pairs = zip(cut_sites[:-1], cut_sites[1:])
        for i, j in pairs:
            expected.append(test_str[i:j])
    expected.append(test_str[cut_sites[-1]:])
    if circular:
        last = expected[-1]
        expected.remove(last)
        expected[0] = last + expected[0]
    assert [str(f) for f in fragments] == expected

@pytest.mark.parametrize('direction', [1,-1])
@pytest.mark.parametrize('query_test_str,expected,template_circular,query_circular', [
    ('ABC', (1, 6), False, False),
    ('ABCD', (1,), False, False),
    ('EABC', (5,), False, False),
    ('234', tuple(), False, False),
    ('BCD', (2,), False, False),
    ('BCD', (2, 7), True, False),
    ('BCDA', (7, ), True, False),
    ('BCD', (2,), False, True),
])
def test_find_iter(query_test_str, expected, template_circular, query_circular, direction):
    query = DoubleLinkedList(data=query_test_str)
    template_test_str = 'DABCDEABC'
    template = DoubleLinkedList(data=template_test_str)
    template.circular = template_circular
    query.circular = query_circular
    found = list(template.find_iter(query, direction=direction))
    found = sorted(found, key=lambda x: template.index_of(x.start))

    # for f in found:
    #     print(template.index_of(f.start))

    assert len(found) == len(expected)
    for i, f in enumerate(found):
        assert f.start is template[expected[i]]
        j = expected[i] + len(query_test_str) - 1
        if template_circular and j >= len(template):
            j -= len(template)
        assert f.end is template[j]
        assert f.span == (expected[i], j)


@pytest.mark.parametrize('direction', [
    1,
    -1
])
@pytest.mark.parametrize('length,expected', [
    (1, ((0, 5), (6, 4), (11, 3))),
    (2, ((0, 5), (6, 4), (11, 3))),
    (3, ((0, 5), (6, 4), (11, 3))),
    (4, ((0, 5), (6, 4), (11, 3))),
    (5, ((0, 5), (6, 4),)),
    (6, ((0, 5),)),
])
def test_find_iter_min_length(length, expected, direction):
    template_test_str = "ABCDEFABCDEABCD"
    template = DoubleLinkedList(template_test_str)

    found = list(template.find_iter(DoubleLinkedList('ABCDEFZZZ'), min_query_length=length, direction=direction))

    print()
    for i, f in enumerate(found):
        assert f.start is template[expected[i][0]]
        print(template[f.span[0]:f.span[1] + 1])
        assert f.span == (expected[i][0], expected[i][0] + expected[i][1])


def test_find_iter_min_length_reverse():
    template_test_str = "CDEF BCDEF ABCDEF"
    template = DoubleLinkedList(template_test_str)

    found = list(template.find_iter(DoubleLinkedList('ZZZABCDEF'), min_query_length=1, direction=-1))

    assert len(found) == 3

    assert found[0].start is template.get(0)
    assert found[0].end is template.get(3)
    assert found[0].span == (0, 3)

    assert found[1].start is template.get(5)
    assert found[1].end is template.get(9)
    assert found[1].span == (5, 9)

    assert found[2].start is template.get(11)
    assert found[2].end is template.get(16)
    assert found[2].span == (11, 16)


def test_find_iter_protocol():

    template_test_str = "ABCDEFGHIJKL"
    template = DoubleLinkedList(template_test_str)

    found = list(template.find_iter(DoubleLinkedList("BCD"), protocol=lambda x,y: False))
    assert len(found) == 0


def test_find_iter_reverse_query():
    template_test_str = "ABCDEFGHIJKL"
    template = DoubleLinkedList(template_test_str)

    found = list(template.find_iter(DoubleLinkedList("ZZIJKL"[::-1]), min_query_length=2,
                                    direction=(-1, 1)))
    assert len(found) == 1
    assert found[0].span == (8, 11)
    assert found[0].query_span == (3, 0)


def test_find_iter_query_span():
    template_test_str = "ABCDEFGHIJKL"
    template = DoubleLinkedList(template_test_str)

    query1 = DoubleLinkedList("BCDEFZZZZ")
    for match in template.find_iter(query1, min_query_length=3):
        assert match.span == (1, 5)
        assert match.query_span == (0, 4)

    query2 = DoubleLinkedList("ZZZZBCDEF")
    for match in template.find_iter(query2, min_query_length=3):
        assert match.query_span == (4, 8)


@pytest.mark.parametrize('circular', [False, True])
def test_reverse(circular):
    test_str = 'XXXXGHHHXHGG'
    l = DoubleLinkedList(data=test_str)
    if circular:
        l.circularize()
    assert str(l) == test_str
    l.reverse()
    assert str(l) == test_str[::-1]


@pytest.mark.parametrize('circular', [
    True,
    False,
    True,
    False,
])
def test_all_nodes(circular, linked_list, test_str):
    linked_list.circular = circular
    assert len(linked_list.all_nodes()) == len(test_str)


def test_remove():
    test_str = 'XXXXGHHHXHGG'
    l = DoubleLinkedList(data=test_str)
    for i in range(len(test_str)):
        l_copy = copy(l)
        l_copy.remove(i)
        if i == 0:
            assert str(l_copy) == test_str[i + 1:]
        elif i == len(test_str):
            assert str(l_copy) == test_str[:-1]
        else:
            assert str(l_copy) == test_str[:i] + test_str[i + 1:]


def test_zip():
    l1 = DoubleLinkedList('abcd')
    l2 = DoubleLinkedList('agasdf')
    for x1, x2 in zip(l1, l2):
        print('{} {}'.format(x1, x2))


def test_zip_with_node():
    l1 = DoubleLinkedList('abcdefghijk')
    l2 = DoubleLinkedList('agasdf')
    l2.circularize()
    l1.circularize()
    for x1, x2 in zip(l2.get(1).fwd(), l1.get(2).fwd()):
        print('{} {}'.format(x1, x2))


class TestLinkedListMagic(object):

    @pytest.mark.parametrize('i', [
        0, 2, 3, -2, 7, -10, None
    ])
    @pytest.mark.parametrize('j', [
        0, 2, 3, 4, 7, -2, -10, None
    ])
    @pytest.mark.parametrize('circular', [
        False
    ])
    def test_splice_magic(self, i, j, circular, linked_list, test_str):
        linked_list.cyclic = circular
        expected = test_str[i:j]
        if test_str[i:j] == '':
            assert linked_list[i:j].is_empty(), "({},{}) should return None".format(i, j)
        else:
            assert str(linked_list[i:j]) == test_str[i:j]

    @pytest.mark.parametrize('i,j', [
        (5, 2),
        (0, 0),
        (2, 2),
        (-2, -4),
    ])
    def test_splice_magic_circular(self, i, j, linked_list, test_str):
        linked_list.circularize()
        if i == j:
            assert linked_list[i:j].is_empty()
            return
        assert str(linked_list[i:j]) == test_str[i:] + test_str[:j]

    def test_(self, linked_list, test_str):
        i = 3
        j = None
        assert str(linked_list[i:j]) == test_str[i:j]

    def test_contains_magic(self, linked_list, test_str):
        for n in linked_list:
            assert n in linked_list
        linked_list_2 = DoubleLinkedList(data=test_str)
        for n in linked_list_2:
            assert not n in linked_list

    def test_splice_copy_magic(self, linked_list, test_str):
        assert str(linked_list[:]) == test_str
        assert linked_list[:] is not linked_list

    @pytest.mark.parametrize('circular', [True, False])
    def test_reverse_magic(self, circular, linked_list, test_str):
        if circular:
            linked_list.circularize()
        assert str(linked_list) == test_str
        assert str(linked_list[::-1]) == test_str[::-1]

    @pytest.mark.parametrize('circular', [True, False])
    def test_len_magic(self, circular, linked_list, test_str):
        linked_list.circular = circular
        assert len(linked_list) == len(test_str)

    def test_add_magic(self):
        l1 = DoubleLinkedList('abc')
        l2 = DoubleLinkedList('def')
        l3 = l1 + l2
        assert str(l3) == 'abcdef'
        assert str(l1) == 'abc'
        assert str(l2) == 'def'
        assert l1 is not l2
        assert l2 is not l3


