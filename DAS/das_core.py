'''
Project: jdna
File: das_core
Author: Justin
Date: 2/6/17

Description: 

'''

from das_contig import *
from das_blast import *

locations = dict(
    database="database",
    designs="designs",
    templates="templates"
)

design_path = 'designs/pmodkan-ho-pact1-z4-er-vpr.gb'

b = BLAST('db', 'templates', design_path, 'database', 'database/results.out', evalue=10.0)
b.makedbfromdir()
b.runblast()
contig_container = b.parse_results(contig_type=Contig.TYPE_BLAST)
contig_container.fuse_circular_fragments()
contig_container.remove_redundant_contigs(include_contigs_contained_within=True, save_linear_contigs=True)

contig_container.filter_perfect()

# generate_random_primers(open_sequence(design_path)[0], 'primers/primers.fasta', num_primers=2)
p = BLAST('primerdb', 'primers', design_path, 'database', 'database/primerresults.out', evalue=1000.0, word_size=15)
p.makedbfromdir()
p.runblast()
primer_container = p.parse_results(contig_type=Contig.TYPE_PRIMER)
primer_container.dump('alignment_viewer/primer_data.json')
primer_container.filter_perfect()

assembly = AssemblyGraph(primers=primer_container, contigs=contig_container)
# assembly.expand_contigs(primer_container.contigs)

assembly.remove_redundant_contigs(include_contigs_contained_within=False, save_linear_contigs=False)
# assembly.break_apart_long_contigs()

contig_container.dump('alignment_viewer/data.json')

paths = assembly.get_all_assemblies(place_holder_size=10)
print paths[0].summary()
print paths[0].circular_gap()
# a.fill_contig_gaps()

# print a.summary()
j5 = J5Assembly(paths[0])
credentials = None
with open('j5_credentials.json') as handle:
    credentials = json.load(handle)
j5.submit(**credentials)
# Parse assembly

# Update templates