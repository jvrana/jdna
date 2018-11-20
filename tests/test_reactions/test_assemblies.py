import pytest
from jdna import Sequence, Reaction

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
    indices = list(range(0, len(seq), int(len(seq)/num_fragments)))
    for i, j in zip(indices[:-1], indices[1:]):
        sequences.append(
            seq[i:j+overlap]
        )
    sequences[0] = Sequence(p[-20:]) + sequences[0]
    return sequences


def test_cyclic_assembly_report(sequences):
    report = Reaction.homology_report(sequences)

    assert len(report['cyclic_paths']) == 1, "There should be a single cyclic assembly path"
    assert len(report['cyclic_paths'][0]) == len(sequences), "The cyclic assembly should contain all 3 fragments"


def test_linear_assembly_report():
    seq1 = Sequence()
    seq2 = Sequence()
    seq3 = Sequence()


    report = Reaction.homology_report([seq1, seq2, seq3])

    assert len(report['linear_paths']) == 2, "There should be a single linear assembly path"
    assert len(report['linear_paths'][0]) == 3, "The linear assembly should contain all 3 fragments"
    assert report['linear_paths'][0] == [0, 1, 2], "The linear assembly should be in the correct order [0, 1, 2]"


def test_linear_assembly():
    final_sequence = "CGTTTTAAGAGCTTGGTGAGCGCTAGGAGTCACTGCCAGGTATCGTTTGAACACGGCATTAGTCAGGGAAGTCATAACACAGTCCTTTCCCGCAATTTTCTTTTTCTATTACTCTTGGCCTCCTCTAGTACACTCTATATTTTTTTATGCCTCGGTAATGATTTTCATTTTTTTTTTTCCACCTAGCGGATGACTCTTTTTTTTTCTTAGCGATTGGCATTATCACATAATGAATTATACATTATATAAAGTAATGTGATTTCTTCGAAGAATATACTAAAAAATGAGCAGGCAAGATAAACGAAGGCAAAG"

    seq1 = Sequence()
    seq2 = Sequence()
    seq3 = Sequence()

    seqs = Reaction.linear_assembly([seq1, seq2, seq3])

    assert len(seqs) == 2
    assert len(seqs[0]) == len(final_sequence)
    assert str(seqs[0]).upper() == final_sequence.upper()


def test_linear_assembly():
    final_sequence = "CGTTTTAAGAGCTTGGTGAGCGCTAGGAGTCACTGCCAGGTATCGTTTGAACACGGCATTAGTCAGGGAAGTCATAACACAGTCCTTTCCCGCAATTTTCTTTTTCTATTACTCTTGGCCTCCTCTAGTACACTCTATATTTTTTTATGCCTCGGTAATGATTTTCATTTTTTTTTTTCCACCTAGCGGATGACTCTTTTTTTTTCTTAGCGATTGGCATTATCACATAATGAATTATACATTATATAAAGTAATGTGATTTCTTCGAAGAATATACTAAAAAATGAGCAGGCAAGATAAACGAAGGCAAAG"

    seq1 = Sequence()
    seq2 = Sequence()
    seq3 = Sequence()

    seqs = Reaction.linear_assembly([seq1, seq2, seq3])

    assert len(seqs) == 2
    assert len(seqs[0]) == len(final_sequence)
    assert str(seqs[0]).upper() == final_sequence.upper()


def test_cyclic_assembly():
    final_sequence = "CGTTTTAAGAGCTTGGTGAGCGCTAGGAGTCACTGCCAGGTATCGTTTGAACACGGCATTAGTCAGGGAAGTCATAACACAGTCCTTTCCCGCAATTTTCTTTTTCTATTACTCTTGGCCTCCTCTAGTACACTCTATATTTTTTTATGCCTCGGTAATGATTTTCATTTTTTTTTTTCCACCTAGCGGATGACTCTTTTTTTTTCTTAGCGATTGGCATTATCACATAATGAATTATACATTATATAAAGTAATGTGATTTCTTCGAAGAATATACTAAAAAATGAGCAGGCAAGATAAACGAAGGCAAAG"

    seq1 = Sequence()
    seq2 = Sequence()
    seq3 = Sequence()

    report = Reaction.homology_report([seq1, seq2, seq3])
    seqs = Reaction.cyclic_assembly([seq1, seq2, seq3])

    assert len(seqs) == 1
    assert len(seqs[0]) == len(final_sequence)
    assert str(seqs[0]).upper() == final_sequence.upper()