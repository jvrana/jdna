"""
Project: jdna
File: region
Author: Justin
Date: 2/21/17

Description: Basic functionality for defining regions of linear, circularized, or reversed
regions of a sequence.

"""
from jdna.regions.context import Context
from jdna.regions.exceptions import RegionError


def force_same_context(error=False):
    """Wrapper that returns False or raises Error if other Region has a
    different context.

    If error==False, then wrapped function returns False. If error=True,
    then raises a RegionError.
    """

    def context_wrapper(fxn):
        def check_context(*args, **kwargs):
            self = args[0]
            other = args[1]
            if not self.same_context(other):
                if error:
                    raise RegionError(
                        "Cannot compare two regions if they have different sequence contexts."
                    )
                else:
                    return False
            return fxn(*args, **kwargs)

        return check_context

    return context_wrapper


class Region:
    """Classifies an abstract region of a sequence. A region is defined by the
    inclusive "start" and "end" positions in context of an arbitrary sequence
    defined by the start_index and length.

    Regions can be circular or linear. For circular regions, negative indicies and indicies greater
    than the length are allowable and will be converted to appropriate indices.

    Direction of the region can be FORWARD or REVERSE or BOTH. For reversed directions, start
    and end positions should be flipped.

    Alternative start_index can be used (DEFAULT: 1) to handle sequences that start at 0 or 1.

    A new Region can be created either by defining the start and end positions or defining the
    start position and length by Region.create(length, circular) ::

        E.g. Linear Region
            length: 9
            start_index: 1
            start: 2
            end: 5
            length = 5-2+1 = 4
            Context:  |-------|
            C_Index:  1.......9
            Region:    2..4

        E.g. Circular Region
            length: 9
            start_index: 1
            start: 8
            end: 2
            length = 4
            Context:  |-------|
            C_Index:  1.......9
            Region:   .2     8.

        E.g. Calling Region(5, 7, 9, True, direction=REVERSE, start_index=1)
                  |-------|         Context
                  123456789
                  <<<<| |<<        Region Direction
                      s e          Start (s) and end (e)
    """

    START_INDEX = 0
    FORWARD = 1
    REVERSE = -1
    BOTH = 2

    def __init__(self, start, end, context, direction=FORWARD, name=""):
        """Annotates some region from a contextual sequence with either
        circular or linear topologies.

        Making a forward region:
            Region(1, 5, context=Context(length=10, circular=False, start_index=1))

        Making a reverse region:
            Region(5, 1, direction=Region.REVERSE, context=Context(length=10, circular=False, start_index=1))

        Invalid regions:
            Region(1, 5, direction=Region.REVERSE, context=Context(length=10, circular=False, start_index=1))
                > This would imply a circular sequence context since the region starts at 1 goes in reverse through an
                origin to index 5.

            Region(5, 1, direction=Region.FORWARD, context=Context(length=10, circular=False, start_index=1))
                This would imply a circular sequence context since the region starts at 5 goes forward through an
                origin to index 1.

        :param start: start index of region
        :type start: int
        :param end: end index of region
        :type end: int
        :param context: context sequence for this region
        :type context: Context
        :param direction: Region.FORWARD or Region.REVERSE
        :type direction: int
        :param name: optional name of region
        :type name: str
        """

        self.name = name
        self.context = context
        self.__start = None
        self.__end = None
        self.start = start
        self.end = end
        if direction == Region.REVERSE:
            # self.__start, self.__end = self.__end, self.__start # no longer reverse start and end points
            pass
        self.__direction = direction  # 1 or -1
        self._validate_direction()
        self._validate_region()
        self.start_extendable = False
        self.end_extendable = False

    @staticmethod
    def create_from_ends(left_end, right_end, context, direction=FORWARD, name=""):
        """Creates a region from the endpoints."""
        s, e = left_end, right_end
        if direction == Region.REVERSE:
            s, e = e, s
        return Region(s, e, context=context, direction=direction, name=name)

    def _validate_direction(self):
        """Validates that the direction key is understood."""
        if self.direction not in [Region.FORWARD, Region.REVERSE, Region.BOTH]:
            raise RegionError(
                "Direction {} not understood. Direction must be Region.FORWARD = {}, Region.REVERSE = {},\
             or Region.BOTH = {}".format(
                    self.direction, Region.FORWARD, Region.REVERSE, Region.BOTH
                )
            )

    def _validate_region(self):
        """Validates that the start and end regions are within the bounds of
        its context."""
        if (
            not self.context.circular
            and self.start > self.end
            and self.direction == Region.FORWARD
        ):
            raise RegionError("START cannot be greater than END for linear regions.")
        if (
            not self.context.circular
            and self.start < self.end
            and self.direction == Region.REVERSE
        ):
            raise RegionError("START cannot be greater than END for linear regions.")

    @property
    def start_index(self):
        """The minimum index for this region (defined by context).

        Alias of bounds_start
        """
        return self.bounds_start

    @property
    def bounds_start(self):
        """The maximum index for the context of this region."""
        return self.context.start

    @property
    def bounds_end(self):
        """The maximum index for the context of this region."""
        return self.context.end

    @property
    def context_length(self):
        """The length of the context this region is in."""
        return self.context.length

    @property
    def length(self):
        """The length of this region."""
        if self._spans_origin():
            return self.context.span(
                self.rp, self.lp
            )  # if it spans origin, reverse span calculation
        else:
            return self.context.span(self.lp, self.rp)

    @property
    def circular(self):
        return self.context.circular

    @property
    def start(self):
        """Gets the start position of the region. Internally reverses the start
        and end positions if direction is reversed.

        :return:
        """
        return self.__start

    @start.setter
    def start(self, x):
        """Sets the start position of the region. Internally reverses the start
        and end positions if direction is reversed.

        :return:
        """
        self.__start = self.context.translate_pos(x)

    @property
    def end(self):
        """Gets the end position of the region. Internally reverses the start
        and end positions if direction is reversed.

        :return:
        """
        return self.__end

    @end.setter
    def end(self, x):
        """Sets the end position of the region. Internally reverses the start
        and end positions if direction is reversed.

        :return:
        """
        self.__end = self.context.translate_pos(x)

    @property
    def rp(self):
        """The right most point.

        The maximum index in the region.
        """
        return max(self.start, self.end)

    @property
    def lp(self):
        """The left most point.

        The minimum index in the region.
        """
        return min(self.start, self.end)

    @lp.setter
    def lp(self, v):
        """Sets the left most point."""
        if self.start < self.end:
            self.start = v
        elif self.end < self.start:
            self.end = v
        elif self.end == self.start:
            self.start = v
            self.end = v

    @rp.setter
    def rp(self, v):
        """Sets the right most point."""
        if self.start > self.end:
            self.start = v
        elif self.end > self.start:
            self.end = v
        elif self.end == self.start:
            self.start = v
            self.end = v

    @property
    def direction(self):
        """Direction of the region."""
        return self.__direction

    @property
    def left_end(self):
        """Returns left end stop.

        e.g. ::

            |--------|
            ^L       ^R

            --|        |---
              ^R       ^L

        :return: index of left region end stop
        :rtype: int
        """
        if self._spans_origin():
            return self.rp
        else:
            return self.lp

    @property
    def right_end(self):
        """Returns left end stop e.g. ::

            |--------|
            ^L       ^R

            --|        |---
              ^R       ^L
        :return: index of right region end stop
        :rtype: int
        """
        if self._spans_origin():
            return self.lp
        else:
            return self.rp

    @left_end.setter
    def left_end(self, x):
        """Sets the left end stop."""
        if self.left_end == self.start:
            self.start = x
        elif self.left_end == self.end:
            self.end = x

    @right_end.setter
    def right_end(self, x):
        """Sets the right end stop."""
        if self.right_end == self.start:
            self.start = x
        elif self.right_end == self.end:
            self.end = x

    def is_forward(self):
        """Whether this region is pointed 'forward'."""
        d = self.direction
        return d in [Region.FORWARD]

    def is_reverse(self):
        """Whether this region is pointed in 'reverse'."""
        d = self.direction
        return d in [Region.REVERSE]

    def _spans_origin(self):
        """Returns whether region spans origin.

        Agnostic as to self.circular.
        """
        if self.start > self.end and self.direction == Region.FORWARD:
            return True
        if self.end > self.start and self.direction == Region.REVERSE:
            return True
        return False

    def _valid_indices(self, inclusive=True):
        """Returns tuples of valid indices within this region. Handles circular
        regions. ::

            E.g. |---------|
                 123456789
                   |>>>|
                   [(3,7)] if inclusive=True
                   [(4,6)] if inclusive=False

            E.g. |---------|
                 123456789
                   ||
                   [(3,4)] if inclusive=True
                   [] if inclusive=False

            E.g. |-------|
                 123456789
                 >>| |>>>>
                   [(1,3),(5,9)] if inclusive=True
                   [(1,2),(6,9)] if inclusive=False

            E.g. |-------|
                 123456789
                 |       |>
                [(1,1),(9,9)] if inclusive=True
                [] if inclusive=False

        :param inclusive: Whether region in inclusive
        :type inclusive: bool
        :return: list of tuples indicating valid index ranges for this region.
        :rtype: list
        """
        if not self._spans_origin():
            v = [(self.lp, self.rp)]
            if not inclusive:
                v = [(self.lp + 1, self.rp - 1)]
            if v[0][0] > v[0][1]:  # there are no valid indices
                v = []
            return v
        else:
            ranges = []
            left = (self.bounds_start, self.lp)
            right = (self.rp, self.context.end)
            if not inclusive:
                left = (self.context.start, self.lp - 1)
                right = (self.rp + 1, self.context.start)
            if left[0] <= left[1]:
                ranges.append(left)
            if right[0] <= right[1]:
                ranges.append(right)
            return ranges

    def within_region(self, pos, inclusive=True):
        """Returns whether a position is within the region. Handles circular
        regions.

        :param pos: index
        :type pos: int
        :param inclusive: whether to be inclusive (or exclusive)
        :type inclusive: bool
        :return: whether index is within region
        :rtype: bool
        """
        v = self._valid_indices(inclusive=inclusive)
        x = [i[0] <= pos <= i[1] for i in v]
        return any(x)

    def sub_region(self, s, e):
        """Creates a sub region with the same context_length and direction as
        this region.

        :param s: start index
        :type s: int
        :param e: end index
        :type e: int
        :return: creates a sub region of this region
        :rtype: Region
        """
        if self.within_region(s, inclusive=True) and self.within_region(
            e, inclusive=True
        ):
            r = self.copy()
            r.start = s
            r.end = e
            r._validate_region()
            return r
        else:
            raise RegionError(
                "Sub region bounds [{}-{}] outside of Region bounds [{}-{}]".format(
                    s, e, self.context.start, self.context.end
                )
            )

    def same_context(self, other):
        """Compares the context properties (length, bounds_start, bounds_ends)
        between two Regions.

        :param other: The other Region
        :type other: Region
        :return: Whether the other region's context sequence is the same length and has same start_index as this Region
        :rtype: bool
        """
        return self.context == other.context

    # TODO: Why not __copy__?
    def copy(self):
        """Creates another region with identical properties.

        :return: Copied region
        :rtype: Region
        """
        # s, e = self.start, self.end
        # if self.direction is Region.REVERSE:
        #     s, e = self.end, self.start
        return self.__class__(
            self.start, self.end, self.context, direction=self.direction, name=self.name
        )

    # TODO: Why does this copy itself??
    @force_same_context(error=True)
    def get_overlap(self, other):
        """Returns a region representing the overlap with another region e.g.
        ::

                |--------|         self
                      |--------|   other
                      |--| << This is returned

        :param other: other Region
        :type other: Region
        :return: Region if there is an overlap, None if there is no overlap
        :rtype: Region
        """
        if self.end_overlaps_with(other):
            r = self.copy()
            r.__direction = Region.FORWARD
            r.start = other.left_end
            r.end = self.right_end
            return r
        else:
            return None

    @force_same_context(error=True)
    def equivalent_location(self, other):
        return other.start == self.start and other.end == self.end

    @force_same_context(error=True)
    def encompasses(self, other, inclusive=True):
        return self.within_region(
            other.start, inclusive=inclusive
        ) and self.within_region(other.end, inclusive=inclusive)

    @force_same_context(error=True)
    def end_overlaps_with(self, other):
        """Whether this region overlaps the next region it this regions end.
        Other needs some kind of overhang to return True. False if context is
        different. ::

             True
                 self   |------|
                 other      |-------|

             False
                 self         |------|
                 other  |-------|

             False
                 self   |------|
                 other    |----|

             True
                 self   |------|
                 other    |-----|
        :param other: other Region
        :type other: Region
        :return: True if region ends overlaps, False if otherwise.
        :rtype: bool
        """

        return self.within_region(
            other.left_end, inclusive=True
        ) and not self.encompasses(other)

    @force_same_context(error=True)
    def get_gap(self, other):
        """Gets the gap region. Returns None if there is no gap or Regions are
        not consecutive. Always are in the forward direction. ::

            Context:        |----------------------|
            This Region:        |-------|
            Other Region:                    |------|
            Gap:                         |==|


            r1:         ----|              |-----
            r2:                  |----|
            r1.get_gap(r2)   |==|
            r2.get_gap(r1)             |==|

        :param other: other Region
        :type other: Region
        :return: gap as a Region
        :rtype: Region
        """
        if self.consecutive_with(other):
            return None
        if not self.within_region(other.lp) or self is other:
            r = self.copy()
            try:
                r.end = other.left_end - 1
                r.start = self.right_end + 1
            except RegionError:
                return None
            return r
        else:
            return None

    @force_same_context(error=True)
    def get_gap_span(self, other):
        """Returns span of gap. Returns 0 if regions are consecutive, negative
        if regions overlap, positive for gaps.

        e.g. ::

            +        |--------|***|-----| (length of three gap)

            -        |------***
                            ***-----|
                            *** (length of three overhang)

            0        |----||----| (0 for consecutive)

            None     |--------|
                        |-----| (None

        :param other:
        :type other:
        :return: span of gap (0 for consecutive, - for overlap, + for gap)
        :rtype: int
        """
        overlap = self.get_overlap(other)
        gap = self.get_gap(other)
        cons = self.consecutive_with(other)
        if cons:
            return 0
        if overlap is not None:
            return -overlap.length
        if gap is not None:
            return gap.length

    # @force_same_context(error=True)
    # def no_overlap(self, other):
    #     """
    #     Returns True if there is no overlap with the other Region
    #
    #     :param other: Other Region
    #     :type other: Region
    #     :return: if there is no overlap
    #     :rtype: bool
    #     """
    #     return not self.within_region(other.start, inclusive=True) \
    #            and not other.within_region(self.end, inclusive=True)

    @force_same_context(error=True)
    def consecutive_with(self, other, ignore_direction=True):
        """Returns whether the right_end is consecutive with the other region's
        left_end ::

            |-------||----|

        :param other: other Region
        :type other: Region
        :return: True if other region is consecutive with this region
        :rtype: bool
        """

        expected_right_end = None
        expected_left_end = None
        try:
            expected_left_end = self.context.translate_pos(self.right_end + 1)
        except RegionError:
            return False
        try:
            expected_right_end = other.context.translate_pos(other.left_end - 1)
        except RegionError:
            return False
        consecutive = (
            self.right_end == expected_right_end and other.left_end == expected_left_end
        )
        return consecutive

    @force_same_context(error=True)
    def fuse(self, other, inplace=True):
        """Fuses this region with other region. If regions are not consecutive,
        raises RegionError ::

            |-----------| (1)
                         |-----------| (2)
            |------------------------| (1 modified)

        :param other: other Region
        :type other: Region
        :return: this Region (self, modified by extension by other Region)
        :rtype: Region
        """
        if self.consecutive_with(other):
            new_region = self
            if not inplace:
                new_region = self.copy()
            new_region.right_end = other.right_end
            return new_region
        else:
            return None
            # raise RegionError(
            #         "Cannot fuse regions [{}-{}] with [{}-{}].".format(self.start, self.end, other.start, other.end))

    def extend_start(self, x):
        """Extends the start by x amount. Retracts start if negative. Raises
        error if retracts past end.

        e.g. extend start by +4 ::

                s          e
                |----------|
            |<<<|----------|
            |--------------|

        :param x:
        :type x:
        :return:
        :rtype:
        """
        self._extend(x, self.start)

    def extend_end(self, x):
        """Extends the end by x amount. Retracts end if negative. Raises error
        if retracts past start.

        e.g. extend end by +4 ::

            s          e
            |----------|
            |----------|>>>|
            |--------------|

        :param x:
        :type x:
        :return:
        :rtype:
        """
        self._extend(x, self.end)

    def extend_left_end(self, x):
        """Extends (or retracts) the left end by x amount."""
        self._extend(x, self.left_end)

    def extend_right_end(self, x):
        """Extends (or retracts) the right end by x amount."""
        self._extend(x, self.right_end)

    def _extend(self, x, end_pos):
        """Extends or retracts by x amount at end_pos."""
        if x < -self.length + 1:
            raise RegionError(
                "Cannot retract end past region start (start: {}, end: {}, x: {})".format(
                    self.start, self.end, x
                )
            )
        if x > self.context.length - self.length:
            raise RegionError(
                "Cannot extend end around origin and past other end (start: {}, end: {}, "
                "x: {})".format(self.start, self.end, x)
            )
        if end_pos == self.right_end:
            self.right_end += x
        elif end_pos == self.left_end:
            self.left_end -= x
        else:
            raise RegionError(
                "Position at {} not at either end. Cannot extend.".format(end_pos)
            )

    def set_forward(self):
        """Reverses direction of region if region is reverse."""
        if not self.is_forward():
            self.reverse_direction()

    def set_reverse(self):
        """Reverses direction of region if region is forward."""
        if self.is_forward():
            self.reverse_direction()

    def reverse_direction(self):
        """Reverses the direction of this region.

        Reverses the start and end indices.
        """
        self.end, self.start = self.start, self.end
        if self.is_forward():
            self.__direction = Region.REVERSE
        else:
            self.__direction = Region.FORWARD
        return self.direction

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __len__(self):
        return self.length

    def __str__(self):
        return "Region(start={start} end={end}, name={name}, direction={direction}, context={context})".format(
            start=self.start,
            end=self.end,
            direction=self.direction,
            context=self.context,
            name=self.name,
        )

    def __repr__(self):
        return str(self)
