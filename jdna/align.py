from Bio import pairwise2, SeqIO
from Bio.Align.Applications import MafftCommandline
from jdna.viewer import SequenceViewer

class Align(object):

    def __init__(self, sequence):
        self._sequence = sequence

    def pairwise(self, other):
        """
        Perform a pairwise alignment using BioPython.

        :param other: the other sequence
        :type other: Sequence | basestring
        :return: list of alignments (tuples)
        :rtype: list
        """
        alignments = pairwise2.align.globalxx(str(self._sequence).upper(), str(other).upper())
        return alignments

    def print_alignment(self, other, max=1):
        """
        Print an alignment with another sequence as a view object. Output will be similar
        to the following:

        .. code::

            > "Alignment" (170bp)


            0         G-GC---G---G----G-C-------------G-TG-----A-T----T---------T--T---ATGTTCATGGACGCCCGGGT
                      | ||   |   |    | |             | ||     | |    |         |  |   ||||||||||||||||||||
                      GAGCCACGCACGTCCCGGCATATTAACTCCAAGCTGGTTCTACTCGGCTGGGCGGGCGTGATTTTATGTTCATGGACGCCCGGGT

            85        ATCAAGGCAGCGGCTCACGCCTCTCCACGCGG--GACAG--GTGAAC--TATC--C-G-ACTAGG---TATCAA-----AG--AC
                      |||||||||||||||   |  | |    | ||  || ||  | |||   |||   | | |  |||   | || |     ||  |
                      ATCAAGGCAGCGGCT---G--T-T----G-GGCAGA-AGAAG-GAA-AATAT-ATCAGGA--AGGCCGT-TC-AGGTTTAGGGA-

        :param other: the other sequence
        :type other: Sequence | basestring
        :param max: maximum number of alignments to display (default=1)
        :type max: int
        :return: None
        :rtype: None
        """
        alignments = self.pairwise(other)
        for a in alignments[:max]:
            mid = ''
            for x1, x2 in zip(a[0], a[1]):
                if '-' not in [x1, x2]:
                    mid += '|'
                else:
                    mid += ' '
            viewer = SequenceViewer([a[0], mid, a[1]], name="Alignment")
            viewer.print()

    def align_sanger_reads(self, abi_filepaths, format='abi', mafft_exe='/usr/local/bin/mafft', verbose=False):
        seqrecords = []
        for filepath in abi_filepaths:
            seqrecords.append(SeqIO.read(filepath, format=format))
        sequences = [self.from_biopython_seqrecord(seq) for seq in seqrecords]
        self.print_mafft([self] + sequences, mafft_exe, verbose)

    def mafft_align(self, sequences, mafft_exe='/usr/local/bin/mafft', verbose=False):
        return self.mafft([self] + sequences, mafft_exe, verbose)

    def print_mafft_align(self, sequences, mafft_exe, verbose=False):
        self.print_mafft([self] + sequences, mafft_exe, verbose)

    @classmethod
    def mafft(cls, sequences, mafft_exe="/usr/local/bin/mafft", verbose=False):
        import tempfile

        in_file = tempfile.mkstemp()[1]
        cls.to_fasta(sequences, in_file)

        mafft_cline = MafftCommandline(mafft_exe, input=in_file)
        if verbose:
            print(mafft_cline)
        mafft_cline.auto = True
        result = mafft_cline()
        if verbose:
            print(result[1])

        out_file = tempfile.mkstemp()[1]
        with open(out_file, 'w') as f:
            f.write(result[0])
        return cls.read_fasta(out_file)

    @classmethod
    def print_mafft(cls, sequences, mafft_exe="/usr/local/bin/mafft", verbose=False):
        seqs = cls.mafft(sequences, mafft_exe, verbose)
        viewer = SequenceViewer(seqs,
                                sequence_labels=['({})'.format(i) for i in range(len(seqs))],
                                apply_indices=list(range(len(seqs))),
                                name='MAFFTA alignment',
                                )
        viewer.metadata['Sequence Names'] = '\n\t' + '\n\t'.join(
            ['{} - {}'.format(i, s.name) for i, s in enumerate(seqs)])
        viewer.print()