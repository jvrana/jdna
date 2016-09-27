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
    pcr(template, 0, 20, 'gtagggatct', len(template)-20, len(template), 'ggtagcagtcag')[0]

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
    fragments = [frag1, frag2.reverse_complement(), frag3, frag4]
    products = Reaction.cyclic_assembly(fragments)
    expected = ''.join([str(x) for x in [oh1, f1, oh2, f2, oh3, f3, oh4, f4]])
    assert len(products[0].search_all(Sequence(sequence=expected))) == 1

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
    frag3.create_feature('feature', 'type', 0, 10)
    frag4.create_feature('feature', 'type', len(frag4)-10, len(frag4)-1)
    fragments = [frag1, frag2.reverse_complement(), frag3, frag4]
    products = Reaction.cyclic_assembly(fragments)
    expected = ''.join([str(x) for x in [oh1, f1, oh2, f2, oh3, f3, oh4, f4]])
    assert len(products[0].search_all(Sequence(sequence=expected))) == 1
    assert len(products[0].get_features()) == 1