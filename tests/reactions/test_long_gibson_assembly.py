from jdna.core import Sequence, Reaction, Convert
from copy import copy
import json
import glob


def test_cyclic_assembly_report():
    seq1 = Sequence(name='seq1',
                         sequence="CCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCA")
    seq2 = Sequence(name='seq2',
                         sequence="GTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCAT")
    seq3 = Sequence(name='seq3',
                         sequence="AGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTC" + "CCAATGCTTAATCAGTGAGGCACCTATCTCA")


    report = Reaction.homology_report([seq1, seq2, seq3])

    assert len(report['cyclic_paths']) == 1, "There should be a single cyclic assembly path"
    assert len(report['cyclic_paths'][0]) == 3, "The cyclic assembly should contain all 3 fragments"


def test_linear_assembly_report():
    seq1 = Sequence(name='seq1',
                         sequence="CCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCA")
    seq2 = Sequence(name='seq2',
                         sequence="GTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCAT")
    seq3 = Sequence(name='seq3',
                         sequence="AGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTC")


    report = Reaction.homology_report([seq1, seq2, seq3])

    assert len(report['linear_paths']) == 2, "There should be a single linear assembly path"
    assert len(report['linear_paths'][0]) == 3, "The linear assembly should contain all 3 fragments"
    assert report['linear_paths'][0] == [0, 1, 2], "The linear assembly should be in the correct order [0, 1, 2]"


def test_linear_assembly():
    final_sequence = "CGTTTTAAGAGCTTGGTGAGCGCTAGGAGTCACTGCCAGGTATCGTTTGAACACGGCATTAGTCAGGGAAGTCATAACACAGTCCTTTCCCGCAATTTTCTTTTTCTATTACTCTTGGCCTCCTCTAGTACACTCTATATTTTTTTATGCCTCGGTAATGATTTTCATTTTTTTTTTTCCACCTAGCGGATGACTCTTTTTTTTTCTTAGCGATTGGCATTATCACATAATGAATTATACATTATATAAAGTAATGTGATTTCTTCGAAGAATATACTAAAAAATGAGCAGGCAAGATAAACGAAGGCAAAG"

    seq1 = Sequence(sequence=final_sequence[:50])
    seq2 = Sequence(sequence=final_sequence[30:130])
    seq3 = Sequence(sequence=final_sequence[100:])

    seqs = Reaction.linear_assembly([seq1, seq2, seq3])

    assert len(seqs) == 2
    assert len(seqs[0]) == len(final_sequence)
    assert str(seqs[0]).upper() == final_sequence.upper()


def test_linear_assembly():
    final_sequence = "CGTTTTAAGAGCTTGGTGAGCGCTAGGAGTCACTGCCAGGTATCGTTTGAACACGGCATTAGTCAGGGAAGTCATAACACAGTCCTTTCCCGCAATTTTCTTTTTCTATTACTCTTGGCCTCCTCTAGTACACTCTATATTTTTTTATGCCTCGGTAATGATTTTCATTTTTTTTTTTCCACCTAGCGGATGACTCTTTTTTTTTCTTAGCGATTGGCATTATCACATAATGAATTATACATTATATAAAGTAATGTGATTTCTTCGAAGAATATACTAAAAAATGAGCAGGCAAGATAAACGAAGGCAAAG"

    seq1 = Sequence(sequence=final_sequence[:50])
    seq2 = Sequence(sequence=final_sequence[30:130])
    seq3 = Sequence(sequence=final_sequence[100:])

    seqs = Reaction.linear_assembly([seq1, seq2, seq3])

    assert len(seqs) == 2
    assert len(seqs[0]) == len(final_sequence)
    assert str(seqs[0]).upper() == final_sequence.upper()