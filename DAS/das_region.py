'''
Project: jdna
File: region
Author: Justin
Date: 2/21/17

Description: Basic functionality for defining regions of linear, circularized, or reversed
regions of a sequence.

'''

class RegionError(Exception):
    def __init__(self, msg):
        self.msg = msg

class Region(object):
    """
    Classifies a region of a sequence. A region is defined by the inclusive "start" and "end"
    positions in context of an arbitrary sequence defined by the start_index and length.

    Regions can be circular or linear. For circular regions, negative indicies and indicies greater
    than the length are allowable and will be converted to appropriate indices.

    Direction of the region can be FORWARD or REVERSE or BOTH. For reversed directions, start
    and end positions should be flipped.

    Alternative start_index can be used (DEFAULT: 1) to handle sequences that start at 0 or 1.

    A new Region can be created either by defining the start and end positions or defining the
    start position and region_span by Region.create(length, circular)

    E.g. Linear Region
        length: 9
        start_index: 1
        start: 2
        end: 5
        region_span = 5-2+1 = 4
        Context:  |-------|
        C_Index:  1.......9
        Region:    2..4

    E.g. Circular Region
        length: 9
        start_index: 1
        start: 8
        end: 2
        region_span = 4
        Context:  |-------|
        C_Index:  1.......9
        Region:   .2     8.
    """
    START_INDEX = 1
    FORWARD = 1
    REVERSE = -1
    BOTH = 2

    def __init__(self, start, end, length, circular, direction=FORWARD, name=None, start_index=START_INDEX):
        """

        :param start:
        :param end:
        :param length:
        :param circular:
        :param start_index:
        """

        self.__length = length  # not allowed to reset length
        self.__circular = circular  # not allowed to reset length
        self.__bounds_start = start_index  # not allowed to reset length
        self.__start = self.translate_pos(start)
        self.__end = self.translate_pos(end)
        if direction == Region.REVERSE:
            self.__start, self.__end = self.__end, self.__start
        self.__direction = direction  # 1 or -1
        if self.direction not in [Region.FORWARD, Region.REVERSE, Region.BOTH]:
            raise RegionError("Direction {} not understood. Direction must be Region.FORWARD = {}, Region.REVERSE = {},\
             or Region.BOTH = {}".format(self.direction, Region.FORWARD, Region.BACKWARD, Region.BOTH))
        self.name = name
        if not self.circular and self.start > self.end:
            raise RegionError("START cannot be greater than END for linear regions.")

    @staticmethod
    def create(length, circular, direction=FORWARD, name=None, start_index=START_INDEX):
        r = Region(start_index, start_index + 1, length, circular, direction=direction, name=name,
                   start_index=start_index)
        return r

    @property
    def length(self):
        return self.__length

    @property
    def circular(self):
        return self.__circular

    @property
    def bounds_start(self):
        return self.__bounds_start

    @property
    def bounds_end(self):
        return self.length + self.bounds_start - 1

    @property
    def start(self):
        """
        Gets the start position of the region. Internally
        reverses the start and end positions if direction is
        reversed.
        :return:
        """
        return self.__start

    @start.setter
    def start(self, x):
        """
        Sets the start position of the region. Internally
        reverses the start and end positions if direction is
        reversed.
        :return:
        """
        self.__start = self.translate_pos(x)

    @property
    def end(self):
        """
        Gets the end position of the region. Internally
        reverses the start and end positions if direction is
        reversed.
        :return:
        """
        return self.__end

    @end.setter
    def end(self, x):
        """
        Sets the end position of the region. Internally
        reverses the start and end positions if direction is
        reversed.
        :return:
        """
        self.__end = self.translate_pos(x)

    @property
    def direction(self):
        return self.__direction

    def set_region_by_start_and_span(self, x, span):
        if span <= 0:
            raise RegionError("Cannot have a span of less than or equal to zero.")
        self.start = x
        self.end = x + span - 1
        return self

    def set_region_by_delta_start_and_span(self, deltax, span):
        self.set_region_by_start_and_span(deltax + self.bounds_start, span)
        return self

    def translate_pos(self, pos):
        if self.circular:
            cleared = False
            while not cleared:
                cleared = True
                if pos > self.bounds_end:
                    pos = pos - self.length
                    cleared = False
                if pos < self.bounds_start:
                    pos = pos + self.length
                    cleared = False
        else:
            if not self.within_bounds(pos, inclusive=True):
                raise RegionError(
                    "Position {} outside of bounds for linear region [{} {}].".format(pos, self.bounds_start,
                                                                                      self.bounds_end))
        return pos

    @property
    def region_span(self):
        if self.end >= self.start:
            return self.end - self.start + 1
        else:
            if self.circular:
                left = self.bounds_end - self.start + 1
                right = self.end - self.bounds_start + 1
                return left + right
            else:
                raise RegionError("START is greater than END for linear region.")

    def within_region(self, pos, inclusive=True):
        s = self.start
        e = self.end
        if not (self.within_bounds(s, inclusive=True) and self.within_bounds(e, inclusive=True)):
            return False
        if s > e:
            if self.circular:
                if inclusive:
                    return s <= pos <= self.bounds_end or \
                           self.bounds_start <= pos <= e
                else:
                    return s < pos <= self.bounds_end or \
                           self.bounds_start <= pos < e
            else:
                raise RegionError("Cannot have START greater than END for non-circular regions.")
        if inclusive:
            return s <= pos <= e
        else:
            return s < pos < e

    def within_bounds(self, pos, inclusive=True):
        if inclusive:
            return self.bounds_start <= pos <= self.bounds_end
        else:
            return self.bounds_start < pos < self.bounds_end

    def sub_region(self, s, e):
        if self.within_region(s, inclusive=True) and self.within_region(e, inclusive=True):
            r = Region(s, e, self.length, self.circular, direction=self.direction, name=self.name,
                       start_index=self.bounds_start)
            return r
        else:
            raise RegionError(
                "Sub region bounds [{}-{}] outside of Region bounds [{}-{}]".format(s, e, self.bounds_start,
                                                                                    self.bounds_end))

    def same_context(self, other):
        return self.circular == other.circular and \
               self.bounds_start == other.bounds_start and \
               self.bounds_end == other.bounds_end

    def copy(self):
        s, e = self.start, self.end
        if self.direction is Region.REVERSE:
            s, e = self.end, self.start
        return Region(s, e, self.length, self.circular,
                      direction=self.direction, name=self.name, start_index=self.bounds_start)

    def get_overlap(self, other):

        if self.end_overlaps_with(other):
            r = self.copy()
            r.start = other.start
            r.end = self.end

            return r
        else:
            return None

    def end_overlaps_with(self, other):
        """
         Whether this region overlaps the next region it this regions end
             True
                 self   |------|
                 other      |-------|

             False
                 self         |------|
                 other  |-------|

             False
                 self   |------|
                 other    |----|
         :param other:
         :return:
         """

        return self.same_context(other) \
               and self.within_region(other.start, inclusive=True) \
               and not self.within_region(other, inclusive=True)


    def get_gap(self, other):
        if self.consecutive_with(other):
            return None
        if self.no_overlap(other):
            r = self.copy()
            r.start = self.end + 1
            r.end = other.start - 1
            return r
        else:
            return None

    def get_gap_degree(self, other):
        overlap = self.get_overlap(other)
        gap = self.get_gap(other)
        cons = self.consecutive_with(other)
        if cons:
            return 0
        if overlap is not None:
            return -overlap.region_span
        if gap is not None:
            return gap.region_span

    def no_overlap(self, other):
        return self.same_context(other) \
                and not self.within_region(other.start, inclusive=True) \
                and not other.within_region(self.end, inclusive=True)


    def consecutive_with(self, other):
        after_self_pos = None
        before_other_pos = None
        try:
            before_other_pos = other.translate_pos(other.start - 1)
        except RegionError as e:
            return False
        try:
            after_self_pos = self.translate_pos(self.end + 1)
        except RegionError as e:
            return False
        fusable = self.same_context(other) and \
                  other.start == after_self_pos and \
                  self.end == before_other_pos
        return fusable


    def fuse(self, other):
        if self.consecutive_with(other):
            self.end += other.region_span
            return self
        else:
            raise RegionError(
                "Cannot fuse regions [{}-{}] with [{}-{}].".format(self.start, self.end, other.start, other.end))