"""
Methods for converting jdna objects to other formats
"""
from Bio import SeqIO
from jdna.viewer import FASTAViewer
from jdna.interface import ClassInterface, Interface


class IOInterface(ClassInterface):
    def from_seqrecord(self, seqrecord):
        return self._cls(
            str(seqrecord.seq), name=seqrecord.id, description=seqrecord.description
        )

    @staticmethod
    def fasta(sequences, out=None):
        fasta_viewer = FASTAViewer(sequences)
        if out:
            with open(out, "w") as f:
                f.write(str(fasta_viewer))
        return fasta_viewer

    @classmethod
    def write_fasta(cls, sequences, out):
        return cls.fasta(sequences, out)

    def parse(self, inpath, format):
        seq_gen = SeqIO.parse(inpath, format=format)
        return [self.from_seqrecord(seq) for seq in seq_gen]

    def read(self, inpath, format):
        seq = SeqIO.read(inpath, format=format)
        return self.from_seqrecord(seq)

    def read_fasta(self, inpath):
        return self.parse(inpath, format="fasta")

    def instance(self, instance):
        return IOInstanceInterface(instance, self)


class IOInstanceInterface(Interface):
    def fasta(self, *sequences, out=None):
        return self._class_interface.fasta([self._inst] + sequences, out=out)

    def write_fasta(self, out):
        return self.fasta(out)

    def json(self, *args, **kwargs):
        return self._inst.json(*args, **kwargs)
