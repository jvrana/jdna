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
    # seqs[0].reverse_complement()
    # print(list(seqs[0].anneal(seqs[1])))
    # print(list(seqs[1].anneal(seqs[0])))
    bindings = Reaction.anneal_sequences(seqs)
    for b in bindings:
        print()
        print(b.primer.global_id)
        print(b.template.global_id)
        print(b.position)

@pytest.mark.parametrize('cyclic', [True, False])
@pytest.mark.parametrize('overhang', [20, 10, 15])
@pytest.mark.parametrize('num_fragments', [1, 3])
@pytest.mark.parametrize('reverse_first', [False, True])
def test_anneal_sequences(seq, num_fragments, generate_sequences, cyclic, overhang, reverse_first):
    sequences = generate_sequences(num_fragments=num_fragments, cyclic=cyclic, overhang=overhang)
    if reverse_first:
        sequences[0].reverse_complement()
    bindings = list(Reaction.anneal_sequences(sequences))

    print(bindings)
    assert len(bindings) == (num_fragments - int(not cyclic)) * 2

    for b in bindings:
        assert b.position.span[1] - b.position.span[0] + 1 == overhang, "Length of overhang should be {}".format(overhang)

    for b in bindings:
        print()
        print(b.position)
        print(b.position.template_anneal)
        print(b.position.query_anneal)


# TODO: deterministic cycles with first global id in front (rotate to index of first)
# TODO: all paths using all pairs shortest paths
def test_interaction_graph(seq, generate_sequences):
    sequences = generate_sequences(4, cyclic=False)
    sequences[0].reverse_complement()
    seq_dict = {x.global_id: x for x in sequences}
    seq_dict.update({-x.global_id: x.copy().reverse_complement() for x in sequences})
    G = Reaction.interaction_graph(sequences)
    for e in G.edges:
        print(e)


def test_linear_assemblies(seq, generate_sequences):
    sequences = generate_sequences(4, cyclic=False)
    sequences[0].reverse_complement()
    seq_dict = {x.global_id: x for x in sequences}
    seq_dict.update({-x.global_id: x.copy().reverse_complement() for x in sequences})
    infos = Reaction.linear_assemblies(sequences)
    print()
    print(len(infos))
    for info in infos:
        print()
        print(info['path'])
        for o in info['binding_positions']:
            print(o)
            print(o.span)
            print(o.query_span)
            print(o.query_start)
            print(o.query_end)
            print(o.query_start in o.primer)
    # if not reverse_first:
    #     for b in bindings:
    #         assert b.position.span[0] == 0
    #         assert b.position.span[-1] == overhang-1
    # else:
    #     for b in bindings[:]:
    #         assert b.position.span[0] == 0
    #         assert b.position.span[-1] == overhang-1
        # assert bindings