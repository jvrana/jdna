'''
Project: jdna
File: test_blast
Author: Justin
Date: 2/23/17

Description: 

'''
from DAS.das_blast import BLAST
from DAS.das_contig import Contig


def test_reindexed_templates():
    b = BLAST('db', 'data/blast_test/reindexed_template', 'data/blast_test/reindexed_design/test.gb', 'data/blast_results', 'data/blast_results/results.out', evalue=10.0, ungapped='',
              penalty=-100, perc_identity=100)
    b.makedbfromdir()
    b.runblast()
    contig_container = b.parse_results(contig_type=Contig.TYPE_BLAST)
    assert len(contig_container.contigs) > 1

    contigs = contig_container.contigs
    contig_container.fuse_circular_fragments()
    contig_container.remove_redundant_contigs(remove_equivalent=True, remove_within=True,
                                              no_removal_if_different_ends=True)
    for c in contigs:
        print c.query.name, c.query.start, c.query.end, c.query.length, c.subject.start, c.subject.end, c.query.length
    assert len(contig_container.contigs) == 1
    c = contig_container.contigs[0]
    assert c.subject.region_span == contig_container.meta.query_length


def test_blast_same():
    b = BLAST('db', 'data/blast_test/reindexed_template', 'data/blast_test/reindexed_design/test_reindexed.gb', 'data/blast_results', 'data/blast_results/results.out', evalue=10.0, ungapped='',
              penalty=-100, perc_identity=100)
    b.makedbfromdir()
    b.runblast()
    contig_container = b.parse_results(contig_type=Contig.TYPE_BLAST)

    contig_container.fuse_circular_fragments()
    contig_container.remove_redundant_contigs(remove_equivalent=True, remove_within=True, no_removal_if_different_ends=False)
    # contig_container.remove_redundant_contigs(remove_within=True)

    contigs = contig_container.contigs
    for c in contigs:
        print c.query.name, c.query.start, c.query.end, c.query.length, c.subject.start, c.subject.end, c.query.length

    # Fuse the fragments and there should only be one contig
    assert len(contig_container.contigs) == 1


def test_pseudo_blast():
    design_path = 'data/blast_test/designs/pmodkan-ho-pact1-z4-er-vpr.gb'
    p = BLAST('primerdb', 'data/blast_test/primers', design_path, '', 'database/primerresults.out')
    print p.perfect_matches()