'''
Project: jdna
File: das_core
Author: Justin
Date: 2/6/17

Description: 

'''

from das_assembler import *

locations = dict(
    database="database",
    designs="designs",
    templates="templates"
)

# def make_primer_alignments(primer_dir, template_path):
#     # TODO: add new alignments to p['alignments'], these cost much more (fragment_construction_cost)
#
#
#     database = makedbfromdir('primers', os.path.join(locations['database'], 'primerdb'))
#
#     # Run blast
#     r = runblast(database,
#                  template_path,
#                  os.path.join(locations['database'], 'primer_results.out'), word_size=15, evalue=1000)
#
#     p = parse_results(r)
#     with open('alignment_viewer/primer_data.json', 'w') as output_handle:
#         json.dump(p, output_handle)
#     return p
#
# def run_alignment():
#     print 'Running alignment against templates'
#     # Make database
#     database = makedbfromdir(locations['templates'], os.path.join(locations['database'], 'db'))
#
#     # TODO: only if circular, then pseudo-circular
#
#     circular = True
#     design_seq = os.path.join(locations['designs'], 'pmodkan-ho-pact1-z4-er-vpr.gb')
#     seq = open_sequence(design_seq)[0]
#     if circular:
#         prefix, suffix = design_seq.split('.')
#         design_seq = prefix + '_pseudocircular.' + suffix
#         save_sequence(design_seq, seq + seq)
#
#     # Run blast
#     r = runblast(database,
#                  design_seq,
#                  os.path.join(locations['database'], 'results.out'), evalue=10)
#
#     # fuse circular alignments
#
#     # open meta data
#     meta = {}
#     with open('database/db.json', 'rU') as handle:
#         meta = json.load(handle)
#
#     # parse results
#     contig_container = parse_results(r, additional_metadata=meta)
#     p['meta']['query_circular'] = circular
#     p['meta']['query_length'] = len(seq)
#
#     # clean up alignments
#     fuse_circular_fragments(p['alignments'])
#
#
#
#     generate_random_primers(seq, 'primers/primers.fasta')
#     primer_results = make_primer_alignments(seq, design_seq)
#     # generate_pcr_products(p, primer_results)
#
#     # TODO: filter out primers with imperfect binding
#     # TODO: filter out alignments with imperfect binding
#
#     # Save results
#     print 'saving results ({})'.format(len(p['alignments']))
#     with open('alignment_viewer/data.json', 'w') as output_handle:
#         json.dump(p, output_handle)
#
#     return p

b = BLAST('db', 'templates', 'designs/pmodkan-ho-pact1-z4-er-vpr.gb', 'database', 'database/results.out', evalue=10.0)
b.makedbfromdir()
b.runblast()
contig_container = b.parse_results()
contig_container.fuse_circular_fragments()


p = BLAST('primerdb', 'primers', 'designs/pmodkan-ho-pact1-z4-er-vpr.gb', 'database', 'database/primerresults.out', evalue=1000.0, word_size=15)
p.makedbfromdir()
p.runblast()
primer_container = p.parse_results()



contig_assembly(contig_container.contigs, primer_container.contigs)
assembly_graph(contig_container)

contig_container.dump('alignment_viewer/data.json')
primer_container.dump('alignment_viewer/primer_data.json')