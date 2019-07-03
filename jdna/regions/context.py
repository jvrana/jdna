from jdna.regions.exceptions import RegionError


class Context(object):
    """ Abstract sequence. Circular or linear. """

    DEFAULT_START_INDEX = 0

    def __init__(
        self, length, circular, name=None, id=None, start_index=DEFAULT_START_INDEX
    ):
        """
        Context constructor

        :param length: length of context
        :type length: int
        :param circular: topology of context; True for circular contexts
        :type circular: bool
        :param start_index: the starting index offset for this context (usually 0 or 1)
        :type start_index: int
        """
        self.name = name
        self.id = id
        self.__start_index = start_index
        self.__length = length
        self.__circular = circular

    @property
    def circular(self):
        """Whether this context is circular"""
        return self.__circular

    @property
    def length(self):
        """The length of this context"""
        return len(self)

    @property
    def start(self):
        """
        Start index of allowable positions on context sequence

        :return:
        :rtype:
        """
        return self.__start_index

    @property
    def end(self):
        """
        End index of allowable positions on context sequence

        :return:
        :rtype:
        """
        return self.length + self.start - 1

    def span(self, x, y):
        """
        Calculates the inclusive distance between two points given the context sequence from left to right.

        ::

        e.g.
              context:    |------|
              positions:   x  y
              distance:    1234
        e.g.
              context:    |------|
              positions:   y  x
              distance:   56  1234

        :param x: position 1
        :type x: int
        :param y: position 2
        :type y: int
        :return:
        :rtype: int
        """
        x = self.translate_pos(x)
        y = self.translate_pos(y)
        mx = max(x, y)
        mn = min(x, y)
        m = mx - mn - 1
        if x > y:
            if self.circular:
                return self.length - m
            else:
                return None
        return int(m + 2)

    def within_bounds(self, pos, inclusive=True):
        """
        Whether a position is withing the bounds of acceptable indices given the context sequence.

        :param pos: position
        :type pos: int
        :param inclusive: whether to be inclusive (or exclusive)
        :type inclusive: bool
        :return: True (within bound) or False (not within bounds)
        :rtype: bool
        """
        if inclusive:
            return self.start <= pos <= self.end
        else:
            return self.start < pos < self.end

    def translate_pos(self, pos):
        """
        Translates the index to an allowable index on a circular sequence context.
        Throws RegionError if context is linear and pos is outside of bounds.

        ::

            Context:  |-------|
            C_Index:  1.......9..11
            Pos:      11
            TransPos: 2

        :param pos:
        :type pos:
        :return:
        :rtype:
        """
        if self.circular:
            cleared = False
            while not cleared:
                cleared = True
                if pos > self.end:
                    pos = pos - self.length
                    cleared = False
                if pos < self.start:
                    pos = pos + self.length
                    cleared = False
        else:
            if not self.within_bounds(pos, inclusive=True):
                raise RegionError(
                    "Position {} outside of bounds for linear region [{} {}].".format(
                        pos, self.start, self.end
                    )
                )
        return pos

    def __eq__(self, other):
        """Whether another context is functionally equivalent"""
        return (
            self.circular == other.circular
            and self.start == other.start
            and self.end == other.end
            and self.length == other.length
        )

    def __len__(self):
        """The length of the context"""
        return self.__length

    def __str__(self):
        return "Context(length={length}, circular={circular}, start_index={start_index})".format(
            length=self.length, circular=self.circular, start_index=self.start
        )

    def __repr__(self):
        return str(self)
