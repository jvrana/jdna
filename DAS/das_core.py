'''
Project: jdna
File: das_core
Author: Justin
Date: 2/6/17

Description: 

'''

from das_assembly import *
from das_blast import *
from das_utilities import *



locations = dict(
    database="database",
    designs="designs",
    templates="templates"
)

design_path = 'designs/pmodkan-ho-pact1-z4-er-vpr.gb'
# design_path = 'designs/hcas9-vpr(sanitized).gb'
# design_path = 'designs/pins-011-pef1a-hcsy4-t2a-dcas9-crpos0-crpos1.gb'


# Sanitize all filename
sanitize_filenames('templates')
sanitize_filenames('designs')

b = BLAST('db', 'templates', design_path, 'database', 'database/results.out', evalue=10.0, ungapped='', penalty=-100, perc_identity=100)
b.makedbfromdir()
b.runblast()
contig_container = b.parse_results(contig_type=Contig.TYPE_BLAST)
contig_container.fuse_circular_fragments()
contig_container.remove_redundant_contigs(remove_equivalent=True, remove_within=True, no_removal_if_different_ends=True)



# contig_container.filter_perfect()

# primers = open_sequence('primers/primers.fasta')
# print str(primers[0])
# exit()


# generate_random_primers(open_sequence(design_path)[0], 'primers/primers.fasta', num_primers=2)
p = BLAST('primerdb', 'primers', design_path, 'database', 'database/primerresults.out', evalue=1000.0, task='blastn-short') #, word_size=15, perc_identity=100, penalty=-100, ungapped='') #, ungapped='', penalty=-100)
p.makedbfromdir()
p.runblast()
primer_container = p.parse_results(contig_type=Contig.TYPE_PRIMER)
primer_container.filter_perfect_subjects()

# primer_container.remove_redundant_contigs(include_contigs_contained_within=False, save_linear_contigs=False)
primer_container.dump('alignment_viewer/primer_data.json')

print "Making Assembly Graph"
assembly = AssemblyGraph(primers=primer_container, contigs=contig_container)
print "Expanding Contigs"
assembly.expand_contigs(primer_container.contigs)

print "Breaking Contigs"
assembly.break_contigs_for_circularization()
assembly.break_contigs_at_endpoints()
assembly.break_contigs_for_circularization()
assembly.break_contigs_at_endpoints()
# assembly.break_contigs_at_endpoints()
# assembly.break_contigs_at_endpoints()
# assembly.remove_redundant_contigs(include_contigs_contained_within=False, save_linear_contigs=False)
print "Removing Redundant Contigs"
assembly.remove_redundant_contigs(remove_equivalent=True, remove_within=False, no_removal_if_different_ends=True)
print "Sorting Contigs"
assembly.sort_contigs()
assembly.dump('alignment_viewer/data.json')

dump_coral_to_json(design_path, 'plasmid_viewer/PlasmidViewer/js/plasmiddata.json', width=1000)

paths = assembly.get_all_assemblies(place_holder_size=10, save_history=True)

f = assembly.assembly_history

# import pandas as pd
# import numpy as np
# import seaborn as sns
# import pylab as plt
# a = []
# for i, row in enumerate(f):
#     g = zip([i]*len(row), row)
#     a += g
# x, y = zip(*a)
# x = [float(X) for X in x]
# y = [float(Y) for Y in y]
# d = pd.DataFrame({'step': x, 'score': y})
# print d
# sns.regplot(x='step', y='score', data=d)
#
# exit()
print paths[0].summary()


paths[0].dump('alignment_viewer/bestassembly.json')

# a.fill_contig_gaps()
# print a.summary()

j5 = J5Assembly(paths[0])
credentials = None
with open('j5_credentials.json') as handle:
    credentials = json.load(handle)
# j5.submit(**credentials)
j5.all()
j5.decode_all_to('assembly_parameters')
r = j5.submit(**credentials)
if 'error_message' in r:
    print r['error_message']
else:
    print "Successful!"
with open('assembly_parameters/results.zip', 'w') as handle:
    handle.write(J5Assembly.decode64(r['encoded_output_file']))
# Parse assembly


# Update templates