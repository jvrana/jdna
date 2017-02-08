'''
Project: jdna
File: das_core
Author: Justin
Date: 2/6/17

Description: 

'''

from contig import *
from das_utilities import *
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
contig_container = b.parse_results()
contig_container.fuse_circular_fragments()
contig_container.remove_redundant_contigs(include_contigs_contained_within=True, save_linear_contigs=True)

contig_container.filter_perfect()

generate_random_primers(open_sequence(design_path)[0], 'primers/primers.fasta', num_primers=2)
p = BLAST('primerdb', 'primers', design_path, 'database', 'database/primerresults.out', evalue=1000.0, word_size=15)
p.makedbfromdir()
p.runblast()
primer_container = p.parse_results(contig_type='primer')
primer_container.dump('alignment_viewer/primer_data.json')
#
#
# primer_container.filter_perfect()

assembly = Assembly(contigs=contig_container.contigs, meta=contig_container.meta.__dict__)
assembly.expand_contigs(primer_container.contigs)
# assembly.break_apart_long_contigs()
# assembly.remove_redundant_contigs(include_contigs_contained_within=False)

contig_container.dump('alignment_viewer/data.json')



paths = assembly.get_all_assemblies(place_holder_size=10)

d = assembly.contig_dictionary
print d.keys()
for p in paths[:10]:
    contig_path = [d[x] for x in p]
    cost = assembly.total_cost(contig_path)
    print cost, assembly.compute_assembly_cost(contig_path), [(c.q_start, c.q_end) for c in contig_path], [Assembly.fragment_cost(c) for c in contig_path]
contig_container.sort_contigs()
# print len(contig_container.contigs)


# make_assembly_graph(contig_container)