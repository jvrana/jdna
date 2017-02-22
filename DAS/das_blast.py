'''
Project: jdna
File: das_blast
Author: Justin
Date: 2/6/17

Description: 

'''
from das_seqio import *
from das_utilities import *
from das_contig import Contig, ContigContainer, ContigContainerMeta
import re
from copy import copy
import tempfile

class BLAST(object):

    # in_dir: templates
    # name: db
    # our_dir: database
    # query_path: designs/.....

    def __init__(self, name, db_in_dir, query_path, db_out_dir, results_out, output_format=7, output_views="qacc sacc score evalue bitscore\
     length nident gapopen gaps qlen qstart qend slen sstart send sstrand qseq sseq", **additional_params):
        """

        :param name:
        :param db_in_dir: templates_directory
        :param query_path: design_goal
        :param db_out_dir: this will become temporary
        :param results_out: this will also become temporary
        :param output_format: don't ever change this
        :param output_views:
        :param additional_params:
        """
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

    def save_query_info(self):
        self.query_circular = False
        self.query_seq = open_sequence(self.query)[0]
        self.query_length = len(self.query_seq)

    def runblast(self):
        cmd_str = "blastn -db {db} -query {query} -out {out} -outfmt {outfmt}"
        for key in self.params:
            cmd_str += ' -{} {}'.format(key, self.params[key])
        r = run_cmd(cmd_str, db=self.db, query=self.query, out=self.out, outfmt=self.outfmt)
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

    def concate_db_to_fsa(self):
        """
        Concatenates sequences in self.db into a .fsa file while saving important metadata
        :return:
        """
        out = self.db + '.fsa'
        fasta, seqs, metadata = concat_seqs(self.db_in_dir, out, savemeta=False)
        self.db_input_metadata = metadata
        self.seqs = seqs
        return out, seqs, metadata

    def makedbfromdir(self):
        """
        Concatenates sequencing files from the db_in_dir directory and
        makes a database from the resulting fasta file
        :return: output_path to blast database
        """
        out, seqs, metadata = self.concate_db_to_fsa()
        return self.makedb(out)



    def perfect_matches(self, rc=True):
        """
        Pseudo-blast for finding perfect sequence matches (i.e. primers)
        :param rc:
        :return:
        """
        out, seqs, metadata = self.concate_db_to_fsa()
        query_seq = open_sequence(self.query)[0].seq
        query_seq = re.sub('[nN]', '.', query_seq)

        fwd_matches = []
        rev_matches = []
        for seq in seqs:
            seq = seq.seq
            try:
                rc_seq = dna_reverse_complement(str(seq))
            except KeyError as e:
                print e
                continue
            seq = re.sub('[nN]', '.', seq)
            rc_seq = re.sub('[nN]', '.', rc_seq)

            for match in re.finditer(str(seq), str(query_seq)):
                c = Contig.create_default_contig()
                c.query.start = match.start()
                c.query.end = match.end()
                c.strand = 'plus'
                c.subject_acc = ''
            if rc:
                for rev_match in re.finditer(str(rc_seq), str(query_seq)):
                    pass

    # TODO: Handle circular subject and queries more cleanly
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
                    contig_dict['subject_circular'] = contig_dict['circular']
                    contig_dict['query_circular'] = False
                    contig_container.add_contig(**contig_dict)
        meta['query_circular'] = self.query_circular
        meta['query_length'] = self.query_length
        meta['query_filename'] = self.query
        meta['query_seq'] = self.query_seq
        meta['contig_seqs'] = self.seqs
        contig_container.meta = ContigContainerMeta(**meta)
        return contig_container


class Aligner(BLAST):
    """
    Wrapper for running BLAST. Creates temporary databases and parses results
    """

    def __init__(self, name, templates, design, **additional_params):

        db_out_dir = tempfile.mkdtemp()
        results_out = tempfile.mktemp(dir=db_out_dir)

        super(Aligner, self).__init__(
            name,
            templates,
            design,
            db_out_dir,
            results_out,
            **additional_params
        )

    def run(self, contig_type):
        self.makedbfromdir()
        self.runblast()
        self.contig_container = self.parse_results(contig_type=contig_type)
        return self.contig_container


