"""
class representing base pairs and their complements
"""

import random
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import re


class Alphabet(object):
    """
    A dictionary class that retrieves complementary base pairs
    """

    __slots__ = ["_chr", "_comp", "__complementary", "__ambiguous"]

    def __init__(self, characters, complementary_characters, ambiguous_characters=None):
        self.__complementary = dict(
            zip(
                characters.lower() + characters.upper(),
                complementary_characters.lower() + complementary_characters.upper(),
            )
        )
        if ambiguous_characters is None:
            ambiguous_characters = {}
        self.__ambiguous = ambiguous_characters

    @classmethod
    def from_biopython_alphabet(
        cls, biopython_alphabet: IUPAC.Alphabet, ambiguous_characters=None
    ):
        seq = Seq(biopython_alphabet.letters)
        seq_letters = str(seq)
        c_seq_letters = str(seq.complement())
        return cls(seq_letters, c_seq_letters, ambiguous_characters)

    def characters(self):
        return self.__complementary.keys()

    def random(self):
        """
        Return random character
        """
        return random.choice(list(self.__complementary.keys()))

    def compare(self, s1, s2, ignore_case=True):
        """
        Compare two sequences. If the second sequence has ambiguous bases, convert
        second sequence to a regex to compare with the first sequence.

        :param s1:
        :type s1:
        :param s2:
        :type s2:
        :return:
        :rtype:
        """
        pattern = "".join(self.__ambiguous.get(b, b) for b in s2)
        if ignore_case:
            match = re.match(pattern, s1, re.IGNORECASE)
        else:
            match = re.match(pattern, s1)
        if match:
            return True
        return False

    def complement(self, basestring):
        return "".join(self[x] for x in basestring)

    def reverse_complement(self, basestring):
        return self.complement(basestring)[::-1]

    def __getitem__(self, item):
        return self.__complementary[item]

    def __contains__(self, item):
        return item in self.__complementary

    # aliases
    rc = reverse_complement
    c = complement


_iupacdict = {
    "M": "[AC]",
    "R": "[AG]",
    "W": "[AT]",
    "S": "[CG]",
    "Y": "[CT]",
    "K": "[GT]",
    "V": "[ACG]",
    "H": "[ACT]",
    "D": "[AGT]",
    "B": "[CGT]",
    "X": "[ACGT]",
    "N": "[ACGT]",
}

iupacdict = {}
for k, v in _iupacdict.items():
    iupacdict[k.lower()] = v.lower()
    iupacdict[k.upper()] = v.upper()

AmbiguousDNA = Alphabet.from_biopython_alphabet(IUPAC.ambiguous_dna, iupacdict)
UnambiguousDNA = Alphabet.from_biopython_alphabet(IUPAC.unambiguous_dna, {})
