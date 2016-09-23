from nose.tools import *

from hydradna import *



def test_reaction():


    def test_pcr(template, fi, fj, foh, ri, rj, roh):
        f = Sequence(sequence=foh + str(template)[fi:fj])
        r = Sequence(sequence=str(template)[ri:rj]+roh).reverse_complement()
        expected = foh + str(template)[fi:rj] + roh
        products = Reaction.pcr(template, f, r)
        assert_equal(str(expected).lower(), str(products[0]).lower())
        return products
    template = Sequence(sequence='AGTCGGCGGATCTATGCTGACTGATGTGTGATGTATGCTTGTGTAGTCGTTGAGTCTGATCTG')
    template.create_feature('test', 'test', 0,3)

    test_pcr(template, 0, 20, 'gtagggatct', 40, len(template), 'ggtagcagtcag')[0].find_feature('test')

    test_pcr(template, 3, 25, 'gtagggatct', len(template)-23, len(template)-3, 'ggtagcagtcag')[0].find_feature('test')

    # test multiple products
    template = Sequence(sequence='AGTCGGCGGATCTATGCTGACTGATGTGTGATGTATGCTTGTGTAGTCGTTGAGTCTGATCTGAGTCGGCGGATCTATGCTGACTGATGTGTGATGTATGCTTGTGTAGTCGTTGAGTCTGATCTG')
    f = Sequence(sequence=str(template)[0:20])
    r = Sequence(sequence=str(deepcopy(template).reverse_complement())[0:20])
    products = Reaction.pcr(template, f, r)
    assert_equal(len(products), 3)

    # test pcr analysis
    # Reaction.pcr_analysis(template, f, r)

    # Test Gibson
    frag1 = Sequence(sequence='AGTCGGCGGATCTATGCTGACTGATGTGTGATGTATGCTTGTGTAGTCGTTGAGTCTGATCTG')
    frag2 = Sequence(sequence='TAGTCGTTGAGTCTGATCTGgtcgtagcgcgagcgttgtggcggattctatatatgttgcGGGGAGTGTTCGGTGCGGTGTTATAG')
    frag3 = Sequence(sequence='GGGGAGTGTTCGGTGCGGTGTTATAGgtcgtagcgcgagcgatcttcttgtggcggattctatatatgttgcAGTCGGCGGATCTATGCTGA')
    frag4 = Sequence(sequence='GGGGAGTGTTCGGTGCGGTGTTATAGgtcgtagcgcgagcgatcttcttgtggcggattctatatatgttgcAGTCGGCGGATCTATGCTGA')

    # Reaction.assembly([frag1, frag2, frag3, frag4])
    # Don't forget to deepcopy all the fragments
    fragments = [frag1, frag2, frag3]
    f1 = Feature('f1', 'type')
    frag1.add_feature(0, 30, f1)

    f2 = Feature('f1', 'type')
    f2.end = 19
    frag3.add_feature(len(frag3)-10, len(frag3)-1, f2)
    fused = Reaction.cyclic_assembly(fragments)
    print fused[0].get_feature_ranges()


    # perumations of all fragments

    # find all cycles

    # Test pcr with features

    # Test pcr where feature is cleaved

    # Test pcr where feature is contained

    # Test insertion with fused feature

    # Test insertion otherwise

    # Test gibson

    # Test gibson features fused