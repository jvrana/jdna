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

# Parts
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

# Sequences
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

# Sequences ZIP

# Target Part Order

# Eugene
pass

# Master Oligo

# Master Plasmid

# Master direct (aka fragments)

# j5 Parameters
import base64


class J5INTERFACE(object):
    def __init__(self):
        self.master_plasmids = []
        self.design_goal = None
        self.master_primers = []
        self.master_fragments = []
        self.eugene = None
        self.parts = None
        self.order = None
        self.parameters = None

    def from_assembly(self, assembly):
        pass

    def parts_from_contigs(self, contigs):
        rows = [
            'Part Source (Sequence Display ID),Part Name,Reverse Compliment?,Start (bp),End (bp),Five Prime Internal Preferred Overhangs?,Three Prime Internal Preferred Overhangs?']
        for c in contigs:
            values = [p.contig_id, coral.seqio.read_dna(p.filename).name, False, p.s_start, p.s_end, '', '']
            values = [str(x) for x in values]
            rows.append(','.join(values))
        self.parts = '\n'.join(rows)
        return self.parts

    def seqs_from_contigs(self, contigs):
        rows = ['Sequence File Name,Format']
        for c in contigs:
            format = "Genbank"
            prefix, suffix = f.split('.')
            if suffix == 'fsa' or suffix == 'fasta':
                format = "Fasta"
            rows.append(','.join(c.filename, format))
        self.sequences = '\n'.join(rows)
        return self.sequences


j5 = J5INTERFACE()
j5.from_assembly(a)
# def from_assembly(self, assembly):
#     # Create encoded parts.csv
#     for c in assembly.contigs:
#         self.add_part(p.id, coral.seqio.read_dna(p.filename), False, p.s_start, p.s_end, '', '')
#         self.add_sequence(p.filename)
#
# def get_target_part_order(self, assembly):
#     for c in assembly.contigs:
#         entry = {
#             '(>Bin) or Part Name': '{}__{}({}-{})'.format(c.contig_id, c.subject_acc, c.s_start, c.s_end),
#             'Direction': 'forward',
#             'Forced Assembly Strategy?': '',
#             'Forced Relative Overhang Position?': '',
#             'Direct Synthesis Firewall?': '',
#             "Extra 5' CPEC overlap bps": '',
#             "Extra 3' CPEC overlap bps": ''
#         }
#         self.order.append(entry)
#     return self.order
#
# def get_zipped_sequences(self):
#     pass
#
# def add_sequence(self, filename):
#     format = "Genbank"
#     prefix, suffix = f.split('.')
#     if suffix == 'fsa' or suffix == 'fasta':
#         format = "Fasta"
#     sequence_file = {
#         "Sequence File": f,
#         "Format": format
#     }
#     self.sequences = sequence_file
#
# def add_part(self, name, source, rc, start, end, five_int, three_int):
#     rc = str(rc).upper()
#     name, source, rc, start, end, five_int, three_int
#     part = part = {
#         "Part Name": name,
#         "Part Source (Sequence Display ID)": source,
#         "Reverse Complement?": rc,
#         "Start (bp)": start,
#         "End (bp)": end,
#         "Five Prime Internal Preferred Overhangs?": five_int,
#         "Three Prime Internal Preferred Overhangs?": three_int
#      }
#     self.parts.append(part)
#     return part
