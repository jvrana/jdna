from Bio import pairwise2, SeqIO
from Bio.Align.Applications import MafftCommandline
from jdna.viewer import SequenceViewer

from jdna.interface import ClassInterface, Interface
import tempfile


class AlignInterface(ClassInterface):
    def mafft(self, sequences, mafft_exe="/usr/local/bin/mafft", verbose=False):
        seqs = self._mafft_sequences(sequences, mafft_exe, verbose)
        viewer = self._mafft_viewer(seqs)
        return viewer

    def _mafft_sequences(
        self, sequences, mafft_exe="/usr/local/bin/mafft", verbose=False
    ):
        in_file = tempfile.mkstemp()[1]
        self._cls.IO.fasta(sequences, in_file)

        mafft_cline = MafftCommandline(mafft_exe, input=in_file)
        if verbose:
            print(mafft_cline)
        mafft_cline.auto = True
        mafft_cline.adjustdirection = True
        result = mafft_cline()
        if verbose:
            print(result[1])

        out_file = tempfile.mkstemp()[1]
        with open(out_file, "w") as f:
            f.write(result[0])
        return self._cls.IO.read_fasta(out_file)

    @staticmethod
    def _mafft_viewer(seqs):
        viewer = SequenceViewer(
            seqs,
            sequence_labels=["({})".format(i) for i in range(len(seqs))],
            apply_indices=list(range(len(seqs))),
            name="MAFFTA alignment",
        )
        viewer.metadata["Sequence Names"] = "\n\t" + "\n\t".join(
            ["{} - {}".format(i, s.name) for i, s in enumerate(seqs)]
        )
        return viewer

    def instance(self, instance):
        return AlignInstanceInterface(instance, self)


class AlignInstanceInterface(Interface):
    def pairwise(self, other):
        """
        Perform a pairwise alignment using BioPython.

        :param other: the other sequence
        :type other: Sequence | basestring
        :return: list of alignments (tuples)
        :rtype: list
        """
        alignments = pairwise2.align.globalxx(
            str(self._inst.upper(), str(other).upper())
        )
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
            mid = ""
            for x1, x2 in zip(a[0], a[1]):
                if "-" not in [x1, x2]:
                    mid += "|"
                else:
                    mid += " "
            viewer = SequenceViewer([a[0], mid, a[1]], name="Alignment")
            viewer.print()

    def sanger_reads(
        self,
        abi_filepaths,
        format="abi",
        mafft_exe="/usr/local/bin/mafft",
        verbose=False,
    ):
        seqrecords = []
        for filepath in abi_filepaths:
            seqrecords.append(SeqIO.read(filepath, format=format))
        sequences = [self._inst.IO.from_seqrecord(seq) for seq in seqrecords]
        return self._class_interface.mafft([self._inst] + sequences, mafft_exe, verbose)

    def mafft(self, sequences, mafft_exe="/usr/local/bin/mafft", verbose=False):
        return self._class_interface.mafft([self._inst] + sequences, mafft_exe, verbose)
