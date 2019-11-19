import random

import pytest

from jdna.reaction import Reaction
from jdna.sequence import Sequence


@pytest.mark.parametrize("direction", [1, -1])
@pytest.mark.parametrize("closeness", [0.1, 0.5, 0.9])
@pytest.mark.parametrize(
    "template_length,primer_length",
    [(200, 10), (200, 30), (200, 20), (200, 100), (10000, 2000)],
)
def test_anneal_small(benchmark, template_length, primer_length, direction, closeness):

    template = Sequence.random(template_length)
    start = int(closeness * len(template))
    end = start + primer_length
    primer = template[start:end]
    if direction == -1:
        primer.reverse_complement()
    benchmark(template.anneal, primer)
    if primer_length < Sequence.DEFAULTS.MIN_ANNEAL_BASES:
        assert not len(list(template.anneal(primer)))
    else:
        assert len(list(template.anneal(primer))) == 1


@pytest.fixture(scope="function")
def generate_sequences():
    def generate_sequences(num_fragments=3, cyclic=True):
        seq = Sequence.random(10000)
        sequences = []
        indices = list(range(0, len(seq), int(len(seq) / (num_fragments))))
        j = 0
        for i, j in zip(indices[:-1], indices[1:]):
            sequences.append(seq[i : j + random.randint(20, 30)])
        sequences.append(seq[j:])
        if cyclic:
            sequences[-1].fuse_in_place(seq[: random.randint(20, 30)])
        for i, s in enumerate(sequences):
            s.name = str(i)
        return sequences

    return generate_sequences


@pytest.mark.parametrize("length", [10, 100, 1000, 10000])
def test_constructor(benchmark, length):
    benchmark(Sequence.random, length)


@pytest.mark.parametrize("num_fragments", [2, 6])
@pytest.mark.parametrize("cyclic", [False, True])
def test_interaction_graph(benchmark, num_fragments, cyclic, generate_sequences):
    sequences = generate_sequences(num_fragments, cyclic=cyclic)
    fxn = Reaction.interaction_graph
    benchmark(fxn, sequences, max_bases=50)


@pytest.mark.parametrize("num_fragments", [2, 6])
@pytest.mark.parametrize("cyclic", [False, True])
def test_benchmark_assemblies(benchmark, num_fragments, cyclic, generate_sequences):
    sequences = generate_sequences(num_fragments, cyclic=cyclic)
    fxn = Reaction.cyclic_assemblies
    if not cyclic:
        fxn = Reaction.linear_assemblies
    benchmark(fxn, sequences, max_bases=50)
