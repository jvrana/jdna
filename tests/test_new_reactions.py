from jdna import *
import random
from copy import copy


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
    frag3 = create_product(f3, oh3, oh4)
    # TODO: WHY do these produce different results for finding homology in the homology graph???? (because of inversions still anneal)

    Reaction.overlap_extension_pcr([frag1, frag2, frag3], 1, 2)
    print len(frag1), len(frag2), len(frag3)
    print "*"*50
    Reaction.overlap_extension_pcr([frag2, frag1, frag3], 1, 2)
    print len(frag1), len(frag2), len(frag3)