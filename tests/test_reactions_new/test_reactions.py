import pytest

from jdna.reaction import Reaction
from jdna.sequence import Sequence
import networkx as nx


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
def generate_sequences(seq):
    """
    A number of overlapping sequences
    :param test_str:
    :type test_str:
    :return:
    :rtype:
    """
    def generate_sequences(num_fragments=3, overhang=20, cyclic=True):
        sequences = []
        indices = list(range(0, len(seq), int(len(seq)/(num_fragments))))
        j = 0
        for i, j in zip(indices[:-1], indices[1:]):
            sequences.append(seq[i:j+overhang])
        sequences.append(seq[j:])
        if cyclic:
            sequences[-1].fuse_in_place(seq[:overhang])
        for i, s in enumerate(sequences):
            s.name = str(i)
        return sequences
    return generate_sequences


def test_(generate_sequences):
    seqs = generate_sequences(3, cyclic=False)
    bindings = Reaction.anneal_sequences(seqs)
    for b in bindings:
        print()
        print(b.primer.global_id)
        print(b.template.global_id)
        print(b.position)

@pytest.mark.parametrize('cyclic', ['cyclic', False])
@pytest.mark.parametrize('overhang_length', [20, 30, 15, 10])
@pytest.mark.parametrize('num_fragments', [1, 3])
@pytest.mark.parametrize('reverse_first', [False, 'reverse_first'])
def test_anneal_sequences(seq, num_fragments, generate_sequences, cyclic, overhang_length, reverse_first):
    cyclic = cyclic == 'cyclic'
    reverse_first = reverse_first == 'reverse_first'
    sequences = generate_sequences(num_fragments=num_fragments, cyclic=cyclic, overhang=overhang_length)

    if reverse_first:
        sequences[0].reverse_complement()
    bindings = list(Reaction.anneal_sequences(sequences))

    if overhang_length >= Sequence.DEFAULTS.MIN_ANNEAL_BASES:
        assert len(bindings) == (num_fragments - int(not cyclic)) * 2
    else:
        assert len(bindings) == 0

    for b in bindings:
        assert b.position.span[1] - b.position.span[0] + 1 == overhang_length, "Length of overhang should be {}".format(overhang_length)

    for i, b in enumerate(bindings):
        if reverse_first:
            if (cyclic and i in [0, 1]) or (not cyclic and i == 0):
                assert str(b.position.template_anneal.reverse_complement()) in str(seq), \
                'Reverse sequence ({}) {} not in sequence\n{}'.format(i, b.position.template_anneal.reverse_complement(), seq)
        else:
            assert str(b.position.template_anneal) in str(seq), 'Sequence ({}) {} not in template\n{}'.format(i, b.position.template_anneal, seq)


# TODO: deterministic cycles with first global id in front (rotate to index of first)
# TODO: all paths using all pairs shortest paths\
@pytest.mark.parametrize('num_fragments', [4])
@pytest.mark.parametrize('cyclic', [False, 'cyclic'])
def test_interaction_graph(num_fragments, cyclic, seq, generate_sequences):
    cyclic = cyclic == 'cyclic'
    sequences = generate_sequences(num_fragments, cyclic=cyclic)
    sequences[0].reverse_complement()
    seq_dict = {x.global_id: x for x in sequences}
    seq_dict.update({-x.global_id: x.copy().reverse_complement() for x in sequences})
    G = Reaction.interaction_graph(sequences)
    assert len(G) == num_fragments * 2
    if cyclic:
        assert len(G.edges) == num_fragments * 2
    else:
        assert len(G.edges) == num_fragments * 2 - 2
    for e in G.edges:
        print(e)


def test_linear_assemblies(seq, generate_sequences):
    sequences = generate_sequences(4, cyclic=False)
    sequences[0].reverse_complement()
    seq_dict = {x.global_id: x for x in sequences}
    seq_dict.update({-x.global_id: x.copy().reverse_complement() for x in sequences})
    infos = Reaction.linear_assemblies(sequences)
    print(infos)
    # if not reverse_first:
    #     for b in bindings:
    #         assert b.position.span[0] == 0
    #         assert b.position.span[-1] == overhang-1
    # else:
    #     for b in bindings[:]:
    #         assert b.position.span[0] == 0
    #         assert b.position.span[-1] == overhang-1
        # assert bindings