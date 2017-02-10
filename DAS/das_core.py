'''
Project: jdna
File: das_core
Author: Justin
Date: 2/6/17

Description: 

'''

from das_contig import *
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

# generate_random_primers(open_sequence(design_path)[0], 'primers/primers.fasta', num_primers=2)
p = BLAST('primerdb', 'primers', design_path, 'database', 'database/primerresults.out', evalue=1000.0, word_size=15)
p.makedbfromdir()
p.runblast()
primer_container = p.parse_results(contig_type='primer')
primer_container.dump('alignment_viewer/primer_data.json')
#
#
# primer_container.filter_perfect()

assembly = AssemblyGraph(primers=primer_container, contigs=contig_container)
assembly.expand_contigs(primer_container.contigs)

assembly.remove_redundant_contigs(include_contigs_contained_within=False, save_linear_contigs=False)
# assembly.break_apart_long_contigs()
# assembly.remove_redundant_contigs(include_contigs_contained_within=False)

contig_container.dump('alignment_viewer/data.json')



paths = assembly.get_all_assemblies(place_holder_size=10)
print len(paths)
for p in paths[:5]:
    print [x.contig_id for x in p.contigs]
    print p.total_cost()
    p.print_contig_path()

# TODO: fill in gaps
# a.add_gapped_contigs


a = paths[0]
filenames = []
for p in a.contigs:
    part = {
        "Part Name": p.contig_id,
        "Part Source (Sequence Display ID)": coral.seqio.read_dna(p.filename).name,
        "Reverse Complement?": "FALSE",
        "Start (bp)": p.s_start,
        "End (bp)": p.s_end,
        "Five Prime Internal Preferred Overhangs?": '',
        "Three Prime Internal Preferred Overhangs?": ''
    }
    filenames.append(p.filename)
    print part

for f in list(set(filenames)):
    format = "Genbank"
    prefix, suffix = f.split('.')
    if suffix == 'fsa' or suffix == 'fasta':
        format = "Fasta"
    sequence_file = {
        "Sequence File": f,
        "Format": format
    }
    print sequence_file

# print a.total_cost()
# print a.contigs[0].start_label
# print a.contigs[0].end_label
# print a.contigs[1].start_label
# print a.contigs[1].end_label
# print d.keys()
# for p in paths[:10]:
#     contig_path = [d[x] for x in p]
#     cost = Assembly.total_cost(contig_path, assembly.meta.query_length)
#     print cost, Assembly.assembly_costs(contig_path, assembly.meta.query_length), [(c.q_start, c.q_end) for c in contig_path], [Assembly.pcr_cost(c) for c in contig_path]
# contig_container.sort_contigs()



# print len(contig_container.contigs)


# make_assembly_graph(contig_container)