import pytest
from jdna import alphabet


@pytest.mark.parametrize('seq,expected', [
    ['AGTC', 'TCAG'[::-1]],
    ['aaaAAAAaaaaa', 'tttTTTTttttt'[::-1]]
])
def test_dna_reverse_complement(seq, expected):
    assert alphabet.DNA.reverse_complement(seq) == expected
    assert alphabet.DNA.rc(seq) == alphabet.DNA.reverse_complement(seq)


@pytest.mark.parametrize('seq,expected', [
    ['AGTC', 'TCAG'],
    ['aaaAAAAaaaaa', 'tttTTTTttttt'],
    ['AGTCN', 'TCAGN'],
])
def test_dna_complement(seq, expected):
    assert alphabet.DNA.complement(seq) == expected
    assert alphabet.DNA.c(seq) == alphabet.DNA.complement(seq)



