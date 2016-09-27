import pytest
from jdna.core import Sequence
from jdna.core import Reaction


# Anneal

def anneal_test(start, end, template):
    primer_seq = str(template)[start:end]
    primer1 = Sequence(sequence=str(template)[start:end])
    ann = Reaction.anneal_threeprime(template, primer1)
    s = ann[0][0]
    l = ann[0][1]
    assert primer_seq == str(template)[s - l:s]


def test_3prime_anneal():
    # test 3' anneal
    start = 10
    end = 30
    template = Sequence(
        sequence='ACGCGGTATGTCTGTCTATTGATGTGGTTCTGATGTGCAGTCTGTCTATTGATGTGTGCGTCATGTACGTTCGCGGCGTATATGGGTATGTCTGTCTATTGTTGCTGTGCTGATGTGCGTGTCTGTATTATGCGGCGA')
    anneal_test(start, end, template)

def test_3primer_anneal_at_end():
    #test 3' anneal
    template = Sequence(
        sequence='ACGCGGTATGTCTGTCTATTGATGTGGTTCTGATGTGCAGTCTGTCTATTGATGTGTGCGTCATGTACGTTCGCGGCGTATATGGGTATGTCTGTCTATTGTTGCTGTGCTGATGTGCGTGTCTGTATTATGCGGCGA')
    start = len(template)-25
    end = len(template)
    anneal_test(start, end, template)


    #test 3' anneal
    start = 0
    end = 20

    template = Sequence(sequence= 'ACGCGGTATGTCTGTCTATTGATGTGGTTCTGATGTGCAGTCTGTCTATTGATGTGTGCGTCATGTACGTTCGCGGCGTATATGGGTATGTCTGTCTATTGTTGCTGTGCTGATGTGCGTGTCTGTATTATGCGGCGA')
    primer_seq = str(template)[start:end]
    primer1 = Sequence(sequence=str(template)[start:end])
    ann = Reaction.anneal_threeprime(template, primer1)
    s = ann[0][0]
    l = ann[0][1]
    assert_true(primer_seq, str(template)[s-l:s])


    #test 3' anneal
    start = 0
    end = 20

    template = Sequence(sequence= 'ACGCGGTATGTCTGTCTATTGATGTGGTTCTGATGTGCAGTCTGTCTATTGATGTGTGCGTCATGTACGTTCGCGGCGTATATGGGTATGTCTGTCTATTGTTGCTGTGCTGATGTGCGTGTCTGTATTATGCGGCGA')
    primer_seq = str(template)[start:end]
    primer1 = Sequence(sequence=str(template)[start:end])
    ann = Reaction.anneal_threeprime(template, primer1, min_bases=30)
    assert_equal(ann, [])

    #test 3' anneal cyclic
    start = 10
    end = 10

    template = Sequence(sequence= 'ACGCGGTATGTCTGTCTATTGATGTGGTTCTGATGTGCAGTCTGTCTATTGATGTGTGCGTCATGTACGTTCGCGGCGTATATGGGTATGTCTGTCTATTGTTGCTGTGCTGATGTGCGTGTCTGTATTATGCGGCGA')
    primer_seq = str(template)[-end:] + str(template)[:start]
    primer1 = Sequence(sequence=primer_seq)
    ann = Reaction.anneal_threeprime(template, primer1)
    assert_equal(ann, [])

    template.make_cyclic()
    ann = Reaction.anneal_threeprime(template, primer1)
    expected = template.slice(len(template)-end, start-1)
    expected = ''.join([str(x) for x in expected])
    assert_equal(expected, primer_seq)

    #test 3' multiplebinding and overhang
    start = 0
    end = 20
    template = Sequence(sequence= 'ACGCGGTATGTCTGTCTATTGATGTGGTTCTGATGTGCAGTCTGTCTATTGATGTGTGCGTCATGTACGTTCGCGGCGTATATGGGTATGTCTGTCTATTGTTGCTGTGCTGATGTGCGTGTCTGTATTATGCGGCGA')
    anneal = 'GTCTGTCTATTGATGTG'
    primer_seq = 'gatagtcgatag' + anneal
    primer1 = Sequence(sequence=primer_seq)
    ann = Reaction.anneal_threeprime(template, primer1)
    assert_true(len(ann) > 1)




    # Anneal

    #test 3' anneal
    start = 10
    end = 30

    template = Sequence(sequence= 'ACGCGGTATGTCTGTCTATTGATGTGGTTCTGATGTGCAGTCTGTCTATTGATGTGTGCGTCATGTACGTTCGCGGCGTATATGGGTATGTCTGTCTATTGTTGCTGTGCTGATGTGCGTGTCTGTATTATGCGGCGA')
    primer_seq = str(template)[start:end]
    primer1 = Sequence(sequence=str(template)[start:end])
    ann = Reaction.anneal_fiveprime(template, primer1)
    s = ann[0][0]
    l = ann[0][1]
    assert_equal(primer_seq, str(template)[s:s+l])


    #test 3' anneal
    start = len(template)-25
    end = len(template)

    template = Sequence(sequence= 'ACGCGGTATGTCTGTCTATTGATGTGGTTCTGATGTGCAGTCTGTCTATTGATGTGTGCGTCATGTACGTTCGCGGCGTATATGGGTATGTCTGTCTATTGTTGCTGTGCTGATGTGCGTGTCTGTATTATGCGGCGA')
    primer_seq = str(template)[start:end]
    primer1 = Sequence(sequence=str(template)[start:end])
    ann = Reaction.anneal_fiveprime(template, primer1)
    s = ann[0][0]
    l = ann[0][1]
    assert_equal(primer_seq, str(template)[s:s+l])


    #test 3' anneal
    start = 0
    end = 20

    template = Sequence(sequence= 'ACGCGGTATGTCTGTCTATTGATGTGGTTCTGATGTGCAGTCTGTCTATTGATGTGTGCGTCATGTACGTTCGCGGCGTATATGGGTATGTCTGTCTATTGTTGCTGTGCTGATGTGCGTGTCTGTATTATGCGGCGA')
    primer_seq = str(template)[start:end]
    primer1 = Sequence(sequence=str(template)[start:end])
    ann = Reaction.anneal_fiveprime(template, primer1)
    s = ann[0][0]
    l = ann[0][1]
    assert_equal(primer_seq, str(template)[s:s+l])


    #test 3' anneal
    start = 0
    end = 20

    template = Sequence(sequence= 'ACGCGGTATGTCTGTCTATTGATGTGGTTCTGATGTGCAGTCTGTCTATTGATGTGTGCGTCATGTACGTTCGCGGCGTATATGGGTATGTCTGTCTATTGTTGCTGTGCTGATGTGCGTGTCTGTATTATGCGGCGA')
    primer_seq = str(template)[start:end]
    primer1 = Sequence(sequence=str(template)[start:end])
    ann = Reaction.anneal_fiveprime(template, primer1, min_bases=30)
    assert_equal(ann, [])

    #test 3' anneal cyclic
    start = 10
    end = 10

    template = Sequence(sequence= 'ACGCGGTATGTCTGTCTATTGATGTGGTTCTGATGTGCAGTCTGTCTATTGATGTGTGCGTCATGTACGTTCGCGGCGTATATGGGTATGTCTGTCTATTGTTGCTGTGCTGATGTGCGTGTCTGTATTATGCGGCGA')
    primer_seq = str(template)[-end:] + str(template)[:start]
    primer1 = Sequence(sequence=primer_seq)
    ann = Reaction.anneal_fiveprime(template, primer1)
    assert_equal(ann, [])

    template.make_cyclic()
    ann = Reaction.anneal_fiveprime(template, primer1)
    expected = template.slice(len(template)-end, start-1, fwd=True)
    expected = ''.join([str(x) for x in expected])
    assert_equal(expected, primer_seq)

    #test 3' multiplebinding and overhang
    start = 0
    end = 20
    template = Sequence(sequence= 'ACGCGGTATGTCTGTCTATTGATGTGGTTCTGATGTGCAGTCTGTCTATTGATGTGTGCGTCATGTACGTTCGCGGCGTATATGGGTATGTCTGTCTATTGTTGCTGTGCTGATGTGCGTGTCTGTATTATGCGGCGA')
    anneal = 'GTCTGTCTATTGATGTG'
    primer_seq = anneal
    primer1 = Sequence(sequence=primer_seq)
    ann = Reaction.anneal_fiveprime(template, primer1)
    assert_true(len(ann) > 1)