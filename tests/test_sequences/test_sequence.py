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
        assert Sequence.random(length).is_empty()
        return
    seq = Sequence.random(length)
    assert len(seq) == length

@pytest.mark.parametrize('length', [
    0,
    1,
    10,
    100,
    2000
])
def test_repr(length):
    seq = Sequence.random(length)
    import textwrap
    print(seq.__repr__())

def test_find_iter_complementary():
    seq = Sequence("CTACAAATTTTACACAGTGGGACGGGCCA")


    matches = list(seq.find_iter(Sequence("TTTT"), direction=-1))
    assert len(matches) == 1

    matches = list(seq.find_iter(Sequence("TTTT"), protocol=lambda x, y: x.complementary(y), direction=-1))
    assert len(matches) == 0

    matches = list(seq.find_iter(Sequence("TGTTTA"), protocol=lambda x, y: x.complementary(y), direction=-1))
    assert len(matches) == 1


@pytest.mark.parametrize('reverse,complement,direction,expected', [
    (False, False, 1, True),
    (False, False, -1, True),
    ('reverse', False, (-1, 1), True),
    ('reverse', False, (1, -1), True),
    ('reverse', False, -1, False),
    ('reverse', False, 1, False),
    (False, 'complement', 1, True),
    (False, 'complement', -1, True),
    (False, 'complement', (1, -1), False),
    (False, 'complement', (-1, 1), False),
    ('reverse', 'complement', 1, False),
    ('reverse', 'complement', -1, False),
    ('reverse', 'complement', (1, -1), True),
    ('reverse', 'complement', (-1, 1), True),
])
@pytest.mark.parametrize('i,j', [
    (10, 20),
    (0, None),
    (10, None)
])
def test_find_iter_reverse_complement(i, j, reverse, complement, direction, expected):
    reverse = reverse == 'reverse'
    complement = complement == 'complement'

    # get slice
    if i is None:
        i = 0
    if j is None:
        j = 30 + i
    template = Sequence('N'*i) + Sequence.random(j-i) + Sequence('N'*10)
    query = template[i:j]

    # manipulate query
    if reverse:
        query.reverse()
    if complement:
        query.complement()

    # get results
    if complement:
        results = list(template.find_iter(query, direction=direction, protocol=lambda x, y: x.complementary(y)))
    else:
        results = list(template.find_iter(query, direction=direction))

    if expected:
        assert len(results) == 1
        res = results[0]
        if i is None:
            expected_i = 0
        else:
            expected_i = i

        if j is None:
            expected_j = len(template) - 1
        else:
            expected_j = expected_i + (j-i) - 1
        assert res.span == (expected_i, expected_j)
    else:
        assert not results


@pytest.mark.parametrize('reverse,complement,expected', [
    # (False, False, True),
    ('reverse', 'complement', True),
])
@pytest.mark.parametrize('i,j', [
    (10, 20),
    (0, None),
    (10, None)
])
def test_anneal_basic(i, j, reverse, complement, expected):
    reverse = reverse == 'reverse'
    complement = complement == 'complement'

    # get slice
    if i is None:
        i = 0
    if j is None:
        j = 30 + i
    template = Sequence('N'*i) + Sequence.random(j-i) + Sequence('N'*10)
    query = template[i:j]

    # manipulate query
    if reverse:
        query.reverse()
    if complement:
        query.complement()

    # results = list(template.anneal(query))
    # assert results
    #
    results = list(template.find_iter(query, direction=(1, -1), min_query_length=10, protocol=lambda x, y: x.complementary(y)))
    print(results)
    results = list(template.anneal_reverse(query))
    print(results)


@pytest.mark.parametrize('reverse_complement', [False, True])
def test_anneal_with_overhang(reverse_complement):

    template = Sequence.random(200)
    anneal = template[20:50]
    overhang = Sequence.random(20) + Sequence('N'*5)
    if reverse_complement:
        anneal.reverse_complement()
    primer = overhang + anneal

    for b in template.anneal(primer):
        print(b.span)
        print(b.query_span)
        assert str(b.five_prime_overhang) == str(overhang)





#
# @pytest.mark.parametrize('reverse_complement', [False, True])
# @pytest.mark.parametrize('overhang', [
#     'A',
#     'ACGTGCTTGCGTGTCGTTGA',
#     ''
# ])
# def test_anneal(overhang, reverse_complement):
#     anneal = "CTACAAATTTACACAGTGGGACGGGCCA"
#     primer = Sequence(overhang + anneal)
#     primer_seq = str(primer)
#     anneal_seq = Sequence(anneal)
#     if reverse_complement:
#         anneal_seq.reverse_complement()
#
#     print(str(anneal_seq))
#     seq = Sequence.random(99) + \
#           Sequence("N") + anneal_seq + Sequence("N") + \
#           Sequence.random(100)
#
#     assert str(primer) == primer_seq
#
#     method = seq.anneal_forward
#     if reverse_complement:
#         method = seq.anneal_reverse
#
#     matches = list(method(primer))
#     # check annealing span
#     assert matches[0].span == (100, 100 + len(anneal)-1)
#
#     if reverse_complement:
#         assert matches[0].query_span == (len(overhang) + len(anneal) - 1, len(overhang))
#     else:
#         assert matches[0].query_span == (len(overhang), len(overhang) + len(anneal) - 1)
#     assert str(matches[0].anneal) == anneal
#
#     # check overhang sequence
#     assert matches[0].three_prime_overhang is None
#     if overhang == '':
#         assert matches[0].five_prime_overhang is None
#     else:
#         assert str(matches[0].five_prime_overhang) == overhang
#
#     # check num matches
#     assert len(matches) == 1
#     assert str(primer) == primer_seq
