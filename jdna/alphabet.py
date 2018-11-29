"""
class representing base pairs and their complements
"""

import random


class Alphabet(object):
    """
    A dictionary class that retrieves complementary base pairs
    """

    __slots__ = ['_chr', '_comp', '__data']

    def __init__(self, characters, complementary_characters):
        self._chr = characters
        self._comp = complementary_characters
        self.__data = dict(zip(
            self._chr.lower() + self._chr.upper(),
            self._comp.lower() + self._comp.upper()
        ))

    def characters(self):
        return self.__data.keys()

    def random(self):
        """
        Return random character
        """
        return random.choice(list(self.__data.keys()))

    def __getitem__(self, item):
        return self.__data[item]

    def __add__(self, other):
        return Alphabet(self._chr + other._chr, self._comp + other._comp)

    def __contains__(self, item):
        return item in self.__data


Ambiguous = Alphabet('N', 'N')
UnambiguousDNA = Alphabet('ATCG', 'TAGC')
AmbiguousDNA = UnambiguousDNA + Ambiguous
DNA = UnambiguousDNA + AmbiguousDNA
UnambiguousRNA = Alphabet('AUCG', 'UAGC')
AmbiguousRNA = UnambiguousRNA + Ambiguous
