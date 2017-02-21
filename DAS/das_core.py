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
from das_j5 import *


locations = dict(
    database="database",
    designs="data/designs",
    templates="data/templates"
)

design_path = os.path.join(locations['designs'], 'pmodkan-ho-pact1-z4-er-vpr.gb')
# design_path = 'designs/hcas9-vpr(sanitized).gb'
# design_path = 'designs/pins-011-pef1a-hcsy4-t2a-dcas9-crpos0-crpos1.gb'


# Sanitize all filename
sanitize_filenames('templates')
sanitize_filenames('designs')

b = BLAST('db', locations['templates'], design_path, 'database', 'database/results.out', evalue=10.0, ungapped='', penalty=-100, perc_identity=100)
b.makedbfromdir()
b.runblast()
contig_container = b.parse_results(contig_type=Contig.TYPE_BLAST)
contig_container.fuse_circular_fragments()
contig_container.remove_redundant_contigs(remove_equivalent=True, remove_within=True, no_removal_if_different_ends=True)



# contig_container.filter_perfect()

# primers = open_sequence("data/primers/primers.fasta')
# print str(primers[0])
# exit()


# generate_random_primers(open_sequence(design_path)[0], "data/primers/primers.fasta', num_primers=2)
p = BLAST('primerdb', "data/primers", design_path, 'database', 'database/primerresults.out', evalue=1000.0, task='blastn-short') #, word_size=15, perc_identity=100, penalty=-100, ungapped='') #, ungapped='', penalty=-100)
p.makedbfromdir()
p.runblast()
primer_container = p.parse_results(contig_type=Contig.TYPE_PRIMER)
primer_container.filter_perfect_subjects()

# primer_container.remove_redundant_contigs(include_contigs_contained_within=False, save_linear_contigs=False)
primer_container.dump('views/alignment_viewer/primer_data.json')


assembly = AssemblyGraph(primers=primer_container, contigs=contig_container)
assembly.expand_contigs(primer_container.contigs)


assembly.break_contigs_at_endpoints()
assembly.break_contigs_at_endpoints()
# assembly.remove_redundant_contigs(include_contigs_contained_within=False, save_linear_contigs=False)
assembly.remove_redundant_contigs(remove_equivalent=True, remove_within=False, no_removal_if_different_ends=True)
assembly.sort_contigs()
assembly.dump('views/alignment_viewer/data.json')






dump_coral_to_json(design_path, 'views/plasmid_viewer/PlasmidViewer/js/plasmiddata.json', width=1000)

paths = assembly.get_all_assemblies(place_holder_size=10, save_history=True)

f = assembly.assembly_history

print paths[0].summary()


paths[0].dump('views/alignment_viewer/bestassembly.json')

# a.fill_contig_gaps()
# print a.summary()




j5 = J5Assembly(paths[0])
contigs = j5.contigs
pairs = zip(contigs[0:-1], contigs[1:])
pairs.append((contigs[-1], contigs[0]))
for l, r in pairs:
    _l = Region(l.query.start, l.query.end, l.query.length/2.0, True,
                direction=Region.FORWARD,
                name=l.query.name,
                start_index=pairs[0][0].query.start)
    _r = Region(r.query.start, r.query.end, r.query.length / 2.0, True,
                direction=Region.FORWARD,
                name=r.query.name,
                start_index=pairs[0][0].query.start)
    print _l.length, _r.length
    print 'span', l.query.region_span, _l.region_span, r.query.region_span, _r.region_span
    print 'regions', _l.start, _l.end, _r.start, _r.end
    print 'queries', l.query.start, l.query.end, r.query.start, r.query.end
    print 'bounds', _l.bounds_start, _l.bounds_end, _r.bounds_start, _r.bounds_end
    print 'gap', _l.get_gap_degree(_r)

# exit()
credentials = None
with open('j5_credentials.json') as handle:
    credentials = json.load(handle)
# j5.submit(**credentials)
j5.all()
j5.decode_all_to('assembly_parameters')
r = j5.submit(**credentials)
if 'error_message' in r:
    print r['error_message']
with open('assembly_parameters/results.zip', 'w') as handle:
    handle.write(J5Assembly.decode64(r['encoded_output_file']))
# Parse assembly


# Update templates

# TODO convert direction fragments, into fragments with new primers
print "DONT FORGET TO convert direction fragments, into fragments with new primers"