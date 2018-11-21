import pytest

from jdna.reaction import Reaction
from jdna.sequence import Sequence


@pytest.fixture(scope='function')
def template():
    return Sequence.random(500)


@pytest.mark.parametrize('o1', [0])
@pytest.mark.parametrize('o2', [0])
@pytest.mark.parametrize('l1', [15])
@pytest.mark.parametrize('l2', [15])
def test_pcr(template, l1, l2, o1, o2):
    pos1 = 100
    pos2 = 200
    p1 = Sequence('N'*o1) + template[pos1:pos1+l1]
    p2 = Sequence('N'*o2) + template[pos2-l2:pos2+1].reverse_complement()

    products = Reaction.pcr(template, p1, p2)
    assert len(products) == 1
    expected_len = o1 + o2 + (pos2-pos1)
    assert len(products[0]) == expected_len
    assert str(products[0]) == 'N'*o1 + str(template)[pos1:pos2] + 'N'*o2


@pytest.fixture(scope='function')
def seq():
    return Sequence.random(300)

@pytest.fixture(scope='function')
def sequences(seq):
    """
    A number of overlapping sequences
    :param test_str:
    :type test_str:
    :return:
    :rtype:
    """
    num_fragments = 3
    overlap = 20
    sequences = []
    indices = list(range(0, len(seq), int(len(seq)/(num_fragments))))
    for i, j in zip(indices[:-1], indices[1:]):
        sequences.append(seq[i:j+overlap])
    sequences.append(seq[j:])
    seq[-20:].fuse_in_place(sequences[0])
    for i, s in enumerate(sequences):
        s.name = str(i)

    for s in sequences[::2]:
        s.reverse_complement()
    # sequences[1].reverse_complement()
    return sequences


def test_anneal_sequences(sequences):
    for b in Reaction.anneal_sequences(sequences):
        print(b[1].name + "->" + b[0].name)
        print(b[-1].span)

        if b[-1].span[0] == 0:
            print("anneal left end")
        if b[-1].span[1] == len(b[0]) - 1:
            print("anneal right end")
        # if b[-1].span[0] == 0:
        #         #     print(b[1].name + "->" + b[0].name)
        #         #     print(b[-1].span)