'''
Project: jdna
File: das_blast
Author: Justin
Date: 2/6/17

Description: 

'''
from das_seqio import *
from das_contig import Contig, ContigContainer, ContigContainerMeta
import re
from copy import copy

class BLAST(object):

    # in_dir: templates
    # name: db
    # our_dir: database
    # query_path: designs/.....

    def __init__(self, name, db_in_dir, query_path, db_out_dir, results_out, output_format=7, output_views="qacc sacc score evalue bitscore\
     length nident gapopen gaps qlen qstart qend slen sstart send sstrand qseq sseq", **additional_params):
        self.query = query_path
        self.save_query_info()
        self.pseudocircular()
        self.save_query_as_fsa()
        self.db_in_dir = db_in_dir
        self.db = os.path.join(db_out_dir, name)
        self.original_query = query_path
        self.outfmt = '\"{} {}\"'.format(output_format, output_views)
        self.name = name
        self.db_out_dir = db_out_dir
        self.out = results_out
        self.db_input_path = None
        self.db_input_metadata = {}
        params = additional_params
        self.params = params

    def save_query_as_fsa(self):
        prefix, suffix = self.query.split('.')

        if suffix == 'gb':
            self.query = gb_to_fsa(self.query, prefix + '.fsa')

    def pseudocircular(self):
        circular = seq_is_circular(self.query)
        seq = self.query_seq
        if circular:
            prefix, suffix = self.query.split('.')
            self.query = prefix + '_pseudocircular.' + suffix
            save_sequence(self.query, seq + seq)
            self.query_circular = True
        print "QUERY", self.query

    def save_query_info(self):
        self.query_circular = False
        self.query_seq = open_sequence(self.query)[0]
        self.query_length = len(self.query_seq)

    def runblast(self):
        cmd_str = "blastn -db {db} -query {query} -out {out} -outfmt {outfmt}"
        for key in self.params:
            cmd_str += ' -{} {}'.format(key, self.params[key])
        run_cmd(cmd_str, db=self.db, query=self.query, out=self.out, outfmt=self.outfmt)
        with open(self.out, 'rU') as handle:
            self.results_raw = handle.read()

    def makedb(self, fasta):
        """
        Makes a blast database from a fasta file
        :param fasta: path to input fasta file
        :param title: title to name database (e.g. "db")
        :param database_output: path to output directory (e.g. "/Users/databases/db")
        :return: path to output directory
        """
        cmd_str = "makeblastdb -in {input} -dbtype nucl -title {title} -out {output}"
        run_cmd(cmd_str, input=fasta, title=self.name, output=self.db)
        self.db_input_path = fasta
        return self.db


    def makedbfromdir(self):
        '''
        Concatenates sequencing files from the db_in_dir directory and
        makes a database from the resulting fasta file
        :return: output_path to blast database
        '''
        out = self.db + '.fsa'
        fasta, seqs, metadata = concat_seqs(self.db_in_dir, out, savemeta=False)
        self.db_input_metadata = metadata
        self.seqs = seqs
        return self.makedb(out)


    def parse_results(self, contig_type=None, delimiter=','):
        '''
        This parses a tabulated blast result
        Contig_Container metadata is defined here
        Contig metadata is saved from the metadata from the seqio.concat_seqs method
        :param contig_type:
        :return:
        '''
        if contig_type is None:
            raise Exception("You must define a contig_type!")
        contig_container = ContigContainer()
        results = self.results_raw
        if results.strip() == '':
            return {'contigs': []}
        g = re.search(
            '# (?P<blast_ver>.+)\n# Query: (?P<query>.+)\\n# Database: (?P<database>.+)\n# Fields: (?P<fields>.+)', results)

        meta = g.groupdict()

        matches = re.findall('\n([^#].*)', results)
        meta['fields'] = re.split('\s*{}\s*'.format(delimiter), meta['fields'])
        # clean up fields
        for i, f in enumerate(meta['fields']):
            meta['fields'][i] = f.replace('.', '') \
                .replace(' ', '_') \
                .replace('%', 'perc')
        for m in matches:
            v = m.split('\t')

            # infer datatypes
            for i, val in enumerate(v):
                try:
                    v[i] = float(val)
                except ValueError as e:
                    pass
                try:
                    v[i] = int(val)
                except ValueError as e:
                    pass

            contig_dict = dict(zip(meta['fields'], v))
            contig_dict['contig_type'] = contig_type
            if self.db_input_metadata:
                if contig_dict['subject_acc'] in self.db_input_metadata:
                    contig_dict.update(self.db_input_metadata[contig_dict['subject_acc']])
                    contig_container.add_contig(**contig_dict)
        meta['query_circular'] = self.query_circular
        meta['query_length'] = self.query_length
        meta['query_filename'] = self.query
        meta['query_seq'] = self.query_seq
        meta['contig_seqs'] = self.seqs
        contig_container.meta = ContigContainerMeta(**meta)
        return contig_container
