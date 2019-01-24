import pytest
from jdna import alphabet


@pytest.mark.parametrize('seq,expected', [
    ['AGTC', 'TCAG'[::-1]],
    ['aaaAAAAaaaaa', 'tttTTTTttttt'[::-1]]
])
def test_dna_reverse_complement(seq, expected):
    assert alphabet.AmbiguousDNA.reverse_complement(seq) == expected
    assert alphabet.AmbiguousDNA.rc(seq) == alphabet.AmbiguousDNA.reverse_complement(seq)


@pytest.mark.parametrize('seq,expected', [
    ['AGTC', 'TCAG'],
    ['aaaAAAAaaaaa', 'tttTTTTttttt'],
    ['AGTCN', 'TCAGN'],
])
def test_dna_complement(seq, expected):
    assert alphabet.AmbiguousDNA.complement(seq) == expected
    assert alphabet.AmbiguousDNA.c(seq) == alphabet.AmbiguousDNA.complement(seq)


def test_ambiguous_characters():

    assert alphabet.AmbiguousDNA.compare('AGTCAG', 'NNNNNN')
    assert alphabet.AmbiguousDNA.compare('AGTCAGA', 'AGTCAGW')
    assert alphabet.AmbiguousDNA.compare('AGTCAGA', 'AGTCAGW')
    assert not alphabet.AmbiguousDNA.compare('AGTCAGC', 'AGTCAGW')