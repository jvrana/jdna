import pytest
import random
from jdna.core import Sequence, Reaction, Feature

def test_pcr():
    def pcr(template, fi, fj, foh, ri, rj, roh):
        f = Sequence(sequence=foh + str(template)[fi:fj])
        r = Sequence(sequence=str(template)[ri:rj]+roh).reverse_complement()
        expected = foh + str(template)[fi:rj] + roh
        products = Reaction.pcr(template, f, r)
        assert str(expected).lower() == str(products[0]).lower()
        return products
    template = Sequence(sequence=''.join([random.choice('atgc') for x in range(200)]))
    random_overhang1 = ''.join([random.choice('atcg') for x in range(20)])
    random_overhang2 = ''.join([random.choice('atcg') for x in range(21)])
    pcr(template, 0, 20, random_overhang1, len(template)-20, len(template), random_overhang2)[0]

def test_cyclic_pcr():
    def pcr(template, fi, fj, foh, ri, rj, roh):
        f = Sequence(sequence=foh + str(template)[fi:fj])
        r = Sequence(sequence=str(template)[ri:rj]+roh).reverse_complement()
        expected = foh + str(template)[fi:] + str(template)[:rj] + roh
        products = Reaction.pcr(template, f, r)
        print foh
        print roh
        print str(template)[ri:rj]
        assert str(expected).lower() == str(products[0]).lower()
        return products
    template = Sequence(sequence=''.join([random.choice('atgc') for x in range(200)]))
    template.make_cyclic()
    random_overhang1 = ''.join([random.choice('atcg') for x in range(20)])
    random_overhang2 = ''.join([random.choice('atcg') for x in range(21)])
    pcr(template, len(template)-40, len(template)-20, random_overhang1, 100, 120, random_overhang2)[0]



def test_gibson():

    def ran_seq(l):
        return ''.join([random.choice('agtc') for x in range(l)])

    def create_product(template, oh1, oh2):
        f = Sequence(sequence=oh1 + str(template)[0:20])
        r = Sequence(sequence=str(template)[len(template)-20:len(template)]+oh2).reverse_complement()
        return Reaction.pcr(template, f, r)[0]

    f1 = Sequence(sequence='AGTCGGCGGATCTATGCTGACTGATGTGTGATGTATGCTTGTGTAGTCGTTGAGTCTGATCTG')
    f2 = Sequence(sequence='TAGTCGTTGAGTCTGATCTGgtcgtagcgcgagcgttgtggcggattctatatatgttgcGGGGAGTGTTCGGTGCGGTGTTATAG')
    f3 = Sequence(sequence='GGGGAGTGTTCGGTGCGGTGTTATAGgtcgtagcgcgagcgatcttcttgtggcggattctatatatgttgcAGTCGGCGGATCTATGCTGA')
    f4 = Sequence(sequence='GGGGAGTGTTCGGTGCGGTGTTATAGgtcgtagcgcgagcgatcttcttgtggcggattctatatatgttgcAGTCGGCGGATCTATGCTGA')

    oh1 = ran_seq(20)
    oh2 = ran_seq(20)
    oh3 = ran_seq(20)
    oh4 = ran_seq(20)
    frag1 = create_product(f1, oh1, oh2)
    frag2 = create_product(f2, oh2, oh3)
    frag3 = create_product(f3, oh3, oh4)
    frag4 = create_product(f4, oh4, oh1)
    fragments = [frag1, frag2, frag3, frag4]
    products = Reaction.cyclic_assembly(fragments)
    expected = ''.join([str(x) for x in [oh1, f1, oh2, f2, oh3, f3, oh4, f4]])
    assert len(products[0].search_all(Sequence(sequence=expected))) == 1


def test_gibson_fail():

    def ran_seq(l):
        return ''.join([random.choice('agtc') for x in range(l)])

    def create_product(template, oh1, oh2):
        f = Sequence(sequence=oh1 + str(template)[0:20])
        r = Sequence(sequence=str(template)[len(template)-20:len(template)]+oh2).reverse_complement()
        return Reaction.pcr(template, f, r)[0]

    f1 = Sequence(sequence='AGTCGGCGGATCTATGCTGACTGATGTGTGATGT')
    f2 = Sequence(sequence='TAGTCGTTGAGTCTGATCTGgtcgtagcgcgagcgttgtggcggattctatatatgttgcGGGGAGTGTTCGGTGCGGTGTTATAG')
    f3 = Sequence(sequence='GGGGAGTGTTCGGTGCGGTGTTATAGgtcgtagcgcgagcgatcttcttgtggcggattctatatatgttgcAGTCGGCGGATCTATGCTGA')
    f4 = Sequence(sequence='GGGGAGTGTTCGGTGCGGTGTTATAGgtcgtagcgcgagcgatcttcttgtggcggattctatatatgttgcAGTCGGCGGATCTATGCTGA')

    oh1 = ran_seq(20)
    oh2 = ran_seq(20)
    oh3 = ran_seq(20)
    oh4 = ran_seq(20)
    frag1 = create_product(f1, oh1, ran_seq(20))
    frag2 = create_product(f2, oh2, oh3)
    frag3 = create_product(f3, oh3, oh4)
    frag4 = create_product(f4, oh4, oh1)
    fragments = [frag1, frag2, frag3, frag4]
    products = Reaction.cyclic_assembly(fragments)
    assert len(products) == 0

def test_gibsons_with_inversions():
    def ran_seq(l):
        return ''.join([random.choice('agtc') for x in range(l)])

    def create_product(template, oh1, oh2):
        f = Sequence(sequence=oh1 + str(template)[0:20])
        r = Sequence(sequence=str(template)[len(template) - 20:len(template)] + oh2).reverse_complement()
        return Reaction.pcr(template, f, r)[0]

    f1 = Sequence(sequence='AGTCGGCGGATCTATGCTGACTGATGTGTGATGTATGCTTGTGTAGTCGTTGAGTCTGATCTG')
    f2 = Sequence(sequence='TAGTCGTTGAGTCTGATCTGgtcgtagcgcgagcgttgtggcggattctatatatgttgcGGGGAGTGTTCGGTGCGGTGTTATAG')
    f3 = Sequence(
        sequence='GGGGAGTGTTCGGTGCGGTGTTATAGgtcgtagcgcgagcgatcttcttgtggcggattctatatatgttgcAGTCGGCGGATCTATGCTGA')
    f4 = Sequence(
        sequence='GGGGAGTGTTCGGTGCGGTGTTATAGgtcgtagcgcgagcgatcttcttgtggcggattctatatatgttgcAGTCGGCGGATCTATGCTGA')

    oh1 = ran_seq(20)
    oh2 = ran_seq(20)
    oh3 = ran_seq(20)
    oh4 = ran_seq(20)
    frag1 = create_product(f1, oh1, oh2)
    frag2 = create_product(f2, oh2, oh3)
    frag3 = create_product(f3, oh3, oh4)
    frag4 = create_product(f4, oh4, oh1)
    feature2 = frag2.create_feature('feature', 'type', 0, 10)
    feature1 = frag1.create_feature('feature', 'type', len(frag2)-11, len(frag2)-1)
    fragments = [frag1, frag2.reverse_complement(), frag3, frag4]
    products = Reaction.cyclic_assembly(fragments)
    expected = ''.join([str(x) for x in [oh1, f1, oh2, f2, oh3, f3, oh4, f4]])
    assert len(products[0].search_all(Sequence(sequence=expected))) == 1
    features = products[0].get_features()
    assert len(features) == 1

def test_gibsons_with_inversions():
    def ran_seq(l):
        return ''.join([random.choice('agtc') for x in range(l)])

    def create_product(template, oh1, oh2):
        f = Sequence(sequence=oh1 + str(template)[0:20])
        r = Sequence(sequence=str(template)[len(template) - 20:len(template)] + oh2).reverse_complement()
        return Reaction.pcr(template, f, r)[0]

    f1 = Sequence(sequence='AGTCGGCGGATCTATGCTGACTGATGTGTGATGTATGCTTGTGTAGTCGTTGAGTCTGATCTG')
    f2 = Sequence(sequence='TAGTCGTTGAGTCTGATCTGgtcgtagcgcgagcgttgtggcggattctatatatgttgcGGGGAGTGTTCGGTGCGGTGTTATAG')
    f3 = Sequence(
        sequence='GGGGAGTGTTCGGTGCGGTGTTATAGgtcgtagcgcgagcgatcttcttgtggcggattctatatatgttgcAGTCGGCGGATCTATGCTGA')
    f4 = Sequence(
        sequence='GGGGAGTGTTCGGTGCGGTGTTATAGgtcgtagcgcgagcgatcttcttgtggcggattctatatatgttgcAGTCGGCGGATCTATGCTGA')

    oh1 = ran_seq(20)
    oh2 = ran_seq(20)
    oh3 = ran_seq(20)
    oh4 = ran_seq(20)
    frag1 = create_product(f1, oh1, oh2)
    frag2 = create_product(f2, oh2, oh3)
    frag3 = create_product(f3, oh3, oh4)
    frag4 = create_product(f4, oh4, oh1)
    feature1 = frag3.create_feature('feature1', 'type', 0, 10)
    feature2 = frag4.create_feature('feature2', 'type', len(frag4)-10, len(frag4)-1)
    fragments = [frag1, frag2.reverse_complement(), frag3, frag4]
    products = Reaction.cyclic_assembly(fragments)
    expected = ''.join([str(x) for x in [oh1, f1, oh2, f2, oh3, f3, oh4, f4]])
    assert len(products) == 1
    assert len(products[0].search_all(Sequence(sequence=expected))) == 1
    features = products[0].get_features()
    feature_names = [f.name for f in features]
    assert feature1.name in feature_names
    assert feature2.name in feature_names

def test_gibson_feature_fusion():
    def ran_seq(l):
        return ''.join([random.choice('agtc') for x in range(l)])


    template = Sequence(sequence='aataaaccagccagccggaagggccgagcgcagaagtggtcctgcaactttatccgcctccatccagtctattaattgttgccgggaagctagagtaagtagttcgccagttaatagtttgcgcaacgttgttgccattgctacaggcatcgtggtgtcacgctcgtcgtttggtatggcttcattcagctccggttcccaacgatcaaggcgagttacatgatcccccatgttgtgcaaaaaagcggttagctccttcggtcctccgatcgttgtcagaagtaagttggccgcagtgttatcactcatggttatggcagcactgcataattctcttactgtcatgccatccgtaagatgcttttctgtgactggtgagtactcaaccaagtcattctgagaatagtgtatgcggcgaccgagttgctcttgcccggcgtcaatacgggata')
    template.create_feature('feature', 'type', 0, len(template)-1)
    template.make_cyclic()
    oh1 = ran_seq(20)
    oh2 = ran_seq(20)

    p1 = Sequence(sequence=str(template)[10:30])
    p2 = Sequence(sequence=str(template)[200:220]).reverse_complement()
    p3 = Sequence(sequence=str(template)[180:200])
    p4 = Sequence(sequence=str(template)[30:50]).reverse_complement()

    frag1 = Reaction.pcr(template, p1, p2)[0]
    frag2 = Reaction.pcr(template, p3, p4)[0]
    print frag2
    fragments = [frag1, frag2]
    products = Reaction.cyclic_assembly(fragments)
    expected = str(template)
    print frag1.get_features()
    print frag2.get_features()
    print template.get_features()
    assert len(products[0].search_all(Sequence(sequence=expected))) == 1

def test_overlap_extension_pcr():
    def ran_seq(l):
        return ''.join([random.choice('agtc') for x in range(l)])

    def create_product(template, oh1, oh2):
        f = Sequence(sequence=oh1 + str(template)[0:20])
        r = Sequence(sequence=str(template)[len(template) - 20:len(template)] + oh2).reverse_complement()
        return Reaction.pcr(template, f, r)[0]


    f1 = Sequence(sequence='AGTCGGCGGATCTATGCTGACTGATGTGTGATGT')
    f2 = Sequence(sequence='TAGTCGTTGAGTCTGATCTGgtcgtagcgcgagcgttgtggcggattctatatatgttgcGGGGAGTGTTCGGTGCGGTGTTATAG')
    f3 = Sequence(
        sequence='GGGGAGTGTTCGGTGCGGTGTTATAGgtcgtagcgcgagcgatcttcttgtggcggattctatatatgttgcAGTCGGCGGATCTATGCTGA')
    f4 = Sequence(
        sequence='GGGGAGTGTTCGGTGCGGTGTTATAGgtcgtagcgcgagcgatcttcttgtggcggattctatatatgttgcAGTCGGCGGATCTATGCTGA')

    oh1 = ran_seq(20)
    oh2 = ran_seq(20)
    oh3 = ran_seq(20)
    oh4 = ran_seq(20)
    frag1 = create_product(f1, oh1, oh2)
    frag2 = create_product(f2, oh2, oh3)

    Reaction.overlap_extension_pcr([frag1, frag2], 1, 2)