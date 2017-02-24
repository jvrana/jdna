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
from pyj5 import *

locations = dict(
    database="database",
    designs="data/designs",
    templates="data/templates"
)

# Choose your design
design_path = os.path.join(locations['designs'], 'pmodkan-ho-pact1-z4-er-vpr.gb')
# design_path = 'designs/hcas9-vpr(sanitized).gb'
# design_path = 'designs/pins-011-pef1a-hcsy4-t2a-dcas9-crpos0-crpos1.gb'

# Sanitize all filename
sanitize_filenames('templates')
sanitize_filenames('designs')

# Run blast against template
b = Aligner('db', locations['templates'], design_path, evalue=10.0, ungapped='', penalty=-100, perc_identity=100)
b.makedbfromdir()
b.runblast()
contig_container = b.parse_results(contig_type=Contig.TYPE_BLAST)
contig_container.fuse_circular_fragments()
contig_container.remove_redundant_contigs(remove_equivalent=True, remove_within=True, no_removal_if_different_ends=True)

# Run a primer blast
p = Aligner('primerdb', 'data/primers', design_path, evalue=1000.0, task='blastn-short')
p.makedbfromdir()
p.runblast()
primer_container = p.parse_results(contig_type=Contig.TYPE_PRIMER)
primer_container.filter_perfect_subjects()

# Dump primers to alignment viewer
primer_container.dump('views/alignment_viewer/primer_data.json')

# Expand and sort contigs
assembly = AssemblyGraph(primers=primer_container, contigs=contig_container)
assembly.expand_contigs(primer_container.contigs)
assembly.break_contigs_at_endpoints()
assembly.break_contigs_at_endpoints()
assembly.remove_redundant_contigs(remove_equivalent=True, remove_within=False, no_removal_if_different_ends=True)
assembly.sort_contigs()

# Dump contigs to alignment viewer
assembly.dump('views/alignment_viewer/data.json')

# Dump design to the plasmid viewer
dump_coral_to_json(design_path, 'views/plasmid_viewer/PlasmidViewer/js/plasmiddata.json', width=1000)

# Calculate Assembly Paths
paths = assembly.get_all_assemblies(place_holder_size=10, save_history=True)

# Best assembly
print paths[0].summary()

# Dump best assembly for plasmid viewer
paths[0].dump('views/alignment_viewer/bestassembly.json')

# Get J5 Credentials
credentials = None
with open('j5_credentials.json') as handle:
    credentials = json.load(handle)

# Define J5 Submission
assembly = paths[0]
targets = []
for c in assembly.contigs:
    seq = J5Sequence(c.filename)
    part = J5Part('{}_{}'.format(c.contig_id, c.subject.name), seq, False, c.subject.start, c.subject.end)
    target = J5Target(part)
    targets.append(target)
submission = J5(targets, J5Parameters(), J5.GIBSON)
submission.run(save_files_to_path='results', **credentials)



# j5.submit(**credentials)

print
print "*"*100
print "DONT FORGET TO convert direction fragments, into fragments with new primers"

