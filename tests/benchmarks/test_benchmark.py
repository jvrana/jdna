import pytest
from jdna.sequence import Sequence

@pytest.mark.parametrize('direction', [1, -1])
@pytest.mark.parametrize('closeness', [0.1,0.5,0.9])
@pytest.mark.parametrize('template_length,primer_length', [
    (200, 10),
    (200, 20),
    (200, 100),
    (10000, 2000)
])
def test_anneal_small(benchmark, template_length, primer_length, direction, closeness):

    template = Sequence.random(template_length)
    start = int(closeness*len(template))
    end = start + primer_length
    primer = template[start:end]
    if direction == -1:
        primer.reverse_complement()
    benchmark(template.anneal, primer)
    assert len(list(template.anneal(primer))) == 1