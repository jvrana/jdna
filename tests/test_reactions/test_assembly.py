import random

import pytest

from jdna.reaction import Assembly
from jdna.reaction import Reaction
from jdna.sequence import Sequence


@pytest.fixture(scope="function")
def seq():
    return Sequence.random(300)


@pytest.fixture(scope="function")
def generate_sequences(seq):
    """A number of overlapping sequences.

    :param test_str:
    :type test_str:
    :return:
    :rtype:
    """

    def generate_sequences(num_fragments=3, overhang=20, cyclic=True):
        sequences = []
        indices = list(range(0, len(seq), int(len(seq) / (num_fragments))))
        j = 0
        for i, j in zip(indices[:-1], indices[1:]):
            sequences.append(seq[i : j + random.randint(20, 30)])
        sequences.append(seq[j:])
        if cyclic:
            sequences[-1].fuse_in_place(seq[: random.randint(20, 30)])
        for i, s in enumerate(sequences):
            s.name = "Sequence {}".format(i)
        random.shuffle(sequences)
        return sequences

    return generate_sequences


@pytest.mark.parametrize("cyclic", [False, True])
def test_assembly_init(generate_sequences, cyclic):
    overhangs = [Sequence.random(20) for i in range(4)]
    templates = [Sequence.random(100) for i in range(4)]
    a = Assembly(templates, overhangs, cyclic=cyclic)
    print(a)


@pytest.mark.parametrize(
    "cyclic", [pytest.param(False, id="linear"), pytest.param(True, id="circular")]
)
@pytest.mark.parametrize("try_reverse_complement", [False, True])
@pytest.mark.parametrize("bind_reverse_complement", [False, True])
def test_interaction_graph(
    cyclic, seq, generate_sequences, try_reverse_complement, bind_reverse_complement
):
    num_fragments = 4
    sequences = generate_sequences(num_fragments, cyclic=cyclic)
    if try_reverse_complement == "try_reverse_complement":
        sequences[1].reverse_complement()
    G = Reaction.interaction_graph(
        sequences, bind_reverse_complement=bind_reverse_complement
    )

    edges = list(G.edges)
    print(len(edges))
    if bind_reverse_complement:
        if cyclic:
            assert len(edges) == 4 * 2
        else:
            assert len(edges) == 3 * 2
    else:
        if cyclic:
            assert len(edges) == 4
        else:
            assert len(edges) == 3


@pytest.mark.parametrize(
    "cyclic", [pytest.param(False, id="linear"), pytest.param(True, id="circular")]
)
def test_linear_paths(generate_sequences, cyclic):
    sequences = generate_sequences(4, cyclic=cyclic)
    G = Reaction.interaction_graph(sequences, bind_reverse_complement=True)
    paths = Reaction.linear_paths(G)

    if not cyclic:
        assert len(paths) == 2
        for p in paths:
            assert len(p) == 4
    else:
        assert len(paths) == 0


@pytest.mark.parametrize(
    "cyclic", [pytest.param(False, id="linear"), pytest.param(True, id="circular")]
)
def test_cyclic_paths(generate_sequences, cyclic):
    sequences = generate_sequences(4, cyclic=cyclic)
    sequences[0].reverse_complement()
    G = Reaction.interaction_graph(sequences, bind_reverse_complement=True)
    paths = Reaction.cyclic_paths(G)
    for p in paths:
        print(p)
    for e in G.edges:
        print(e)
    if cyclic:
        assert len(paths) == 2
        for p in paths:
            assert len(p) == 4
    else:
        assert len(paths) == 0


@pytest.mark.parametrize("reverse_complement", [[], [0], [1], [0, 1], [0, 2]])
def test_linear_assemblies(seq, generate_sequences, reverse_complement):
    sequences = generate_sequences(4, cyclic=False)
    for rc in reverse_complement:
        sequences[rc].reverse_complement()
    assemblies = Reaction.linear_assemblies(sequences)

    assert len(assemblies) == 2

    for a in assemblies:
        p = a.product
        assert not p.cyclic
        assert len(p) == len(seq)
        assert str(p) == str(seq) or str(p.copy().reverse_complement()) == str(seq)
        print(a)


@pytest.mark.parametrize("reverse_complement", [[], [0], [1], [2], [0, 1], [0, 2]])
def test_cyclic_assemblies(seq, generate_sequences, reverse_complement):
    sequences = generate_sequences(4, cyclic=True)
    for rc in reverse_complement:
        sequences[rc].reverse_complement()

    assemblies = Reaction.cyclic_assemblies(sequences, max_bases=50)

    assert len(assemblies) == 2

    for s in sequences:
        print(s)
    print()

    for a in assemblies:
        p = a.product
        assert p.cyclic
        expected = seq.copy().circularize()
        rc_expected = expected.copy().reverse_complement()
        g = [p.compare(expected), p.compare(rc_expected)]
        assert any(g)
        a.print(width=50)


def test_cyclic_assemblies_num_fragments(seq, generate_sequences):
    sequences = generate_sequences(10, cyclic=True, overhang=30)
    assemblies = Reaction.cyclic_assemblies(sequences, max_bases=50)
    for a in assemblies:
        print(a)


def test_hard_coded_assembly():
    seqs = [
        Sequence(
            "GTCGGCGGGACCAGGGAGTTTAAACAGGATTGATAATGTAATAGGATCAATGAATATTAACATTAGGTGCTGTGGGTGGCGCTGGAGAAAACCTTCGTATCGGC",
            name="seq1",
        ),
        Sequence(
            "GCCGATACGAAGGTTTTCTCCAGCGAGTTTATCATTATCAGGTTTTGGGACGCTCGAAGGCTTTAATTTGCTTCAATAAAGGAGCGAGCACCCG",
            name="seq2",
        ),
    ]
    Reaction.interaction_report(
        Reaction.interaction_graph(
            seqs, bind_reverse_complement=True, min_bases=10, max_bases=200
        )
    )
    # assemblies = Reaction.cyclic_assemblies(seqs)
    # for a in assemblies:
    #     print(a.product)
    #     a.print(width=100)
