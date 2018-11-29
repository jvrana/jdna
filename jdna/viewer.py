import functools
import itertools
from networkx import nx
from collections import OrderedDict

class StringColumn(object):
    """Class for managing string columns"""

    FILL = ' '

    def __init__(self, strings=None):
        self._strings = []
        self._length = 0
        if strings:
            max_length = max([len(s) for s in strings])
            self._length = max_length
            for s in strings:
                self.append_string(self.right_fill(s))

    @property
    def length(self):
        return self._length

    @property
    def strings(self):
        return self._strings[:]

    def indent(self, num):
        sc_copy = self.copy()
        for i, s in enumerate(sc_copy.strings):
            sc_copy._strings[i] = sc_copy.FILL * num + sc_copy._strings[i]
        sc_copy._length += num
        return sc_copy

    def indent_right(self, num):
        sc_copy = self.copy()
        for i, s in enumerate(sc_copy.strings):
            sc_copy._strings[i] = sc_copy._strings[i] + sc_copy.FILL * num
        sc_copy._length += num
        return sc_copy

    def right_fill(self, string):
        return string + self.FILL * (self.length - len(string))

    def prepend_string(self, new_string):
        if len(new_string) > self.length:
            self._length = len(new_string)
        self._strings.insert(0, self.right_fill(new_string))

    def append_string(self, new_string):
        if len(new_string) > self.length:
            self._length = len(new_string)
        self._strings.append(self.right_fill(new_string))

    def add_prefix(self, prefix):
        for i, s in self.strings:
            self.strings[i] = prefix + self.strings[i]

    def stack(self, *others):
        sc = self[:]
        for other in others:
            for string in other.strings:
                sc.append_string(string)
        return sc

    def __contains__(self, item):
        return any([item in s for s in self.strings])

    def __add__(self, other):
        if isinstance(other, str):
            other = StringColumn([other])
        sc = StringColumn(self.strings)
        diff = len(sc.strings) - len(other.strings)
        if diff > 0:
            for i in range(diff):
                other.prepend_string('')
        elif diff < 0:
            for i in range(-diff):
                sc.prepend_string('')

        new_sc = StringColumn([])
        for this_string, other_string in zip(sc.strings, other.strings):
            new_sc.append_string(this_string + other_string)
        return new_sc

    def copy(self):
        return self.__copy__()

    def strip_indices(self):
        n1 = 0
        n2 = 0
        for x in self[:]:
            if all([_x == ' ' for _x in x]):
                n1 += 1
            else:
                break
        for x in self[::-1]:
            if all([_x == ' ' for _x in x]):
                n2 += 1
            else:
                break
        return n1, len(self) - n2

    def strip(self):
        n1, n2 = self.strip_indices()
        return self[n1:n2]

    def __copy__(self):
        return self.__class__(self.strings)

    def __getitem__(self, key):
        strings = [s.__getitem__(key) for s in self.strings]
        return self.__class__(strings)

    #     def __setitem__(self, key, items):
    #         if not len(items) == len(self.strings):
    #             raise TypeError("Value must have {} items".format(len(self.strings)))
    #         for string, item in zip(self.strings, items):
    #             string[key] = item

    def __eq__(self, other):
        return str(self) == str(other)

    def __iter__(self):
        return zip(*self.strings)

    def __len__(self):
        return self.length

    def __str__(self):
        return '\n'.join(self.strings)

    def __repr__(self):
        return str(self)

    @classmethod
    def condense(cls, rows):
        """
        Condense a list of :class:`StringColumn` into the minimum number of StringColumns comprising of columns stripped
        of white space. Briefly, this is similar to the following procedure:

        .. code-block::

            input = [
                'label         ',
                '       label2 ',
                '      label3  '
            ]

            # >> CONDENSE

            output = [
                'label  label2 ',
                '      label3  '
            ]

        :param rows:
        :type rows:
        :return:
        :rtype:
        """
        segments = []
        indexed_segments = []
        previous_end = 0
        for row in rows:
            start, end = row.strip_indices()
            word = (start, end, row.strip())
            if word not in segments:
                segments.append(tuple(list(word)))
                indexed_segments.append((start, end, row.strip(), previous_end))
                previous_end += 1

        # create a graph of non-overlapping segments
        nonoverlap_graph = nx.Graph()
        for w in indexed_segments:
            nonoverlap_graph.add_node(w[-1])
        for segment1, segment2 in itertools.combinations(indexed_segments, 2):
            start1, end1, _, index1 = segment1
            start2, end2, _, index2 = segment2
            if start1 < start2 or start1 > end2:
                if start2 < start1 or start2 > end1:
                    if end1 < start2:
                        nonoverlap_graph.add_edge(index1, index2)
                    else:
                        nonoverlap_graph.add_edge(index2, index1)

        # find minimum number of cliques that covers the graph (clique covering)
        subgraph = nonoverlap_graph.subgraph(nonoverlap_graph.nodes)
        cliques = []
        while len(subgraph):
            max_clique = list(nx.find_cliques(subgraph))[0]
            cliques.append(max_clique)
            remaining = set(subgraph.nodes).difference(set(max_clique))
            subgraph = nonoverlap_graph.subgraph(list(remaining))

        condensed_rows = []
        for clique in cliques:
            string_column = cls()
            clique_segments = [indexed_segments[s] for s in clique]
            clique_segments = sorted(clique_segments, key=lambda seg: seg[1])

            previous_end = 0
            for segment in clique_segments:
                start, end, seg_str_col, _ = segment
                string_column += seg_str_col.indent(start - previous_end)
                previous_end = end
            condensed_rows.append(string_column)
        return condensed_rows


def chunkify(iterable, n):
    """Break an interable into chunks of size at most 'n'"""
    chunk = None
    for i, x in enumerate(iterable):
        if i % n == 0:
            if chunk is not None:
                yield chunk
            chunk = []
        chunk.append(x)
    yield chunk


def to_lines(string, width):
    """Converts a string to lines of length <= width"""
    lines = []
    for i in range(0, len(string), width):
        lines.append(string[i:i + width])
    return lines


def prepend_lines(lines, label_iterable, indent, fill=' ', align='<'):
    """
    Prepend lines with a label

    :param lines: lines to prepend
    :type lines: list
    :param indent: number of spaces between start of label and start of line
    :type indent: int
    :param fill: default ' '
    :type fill: what to fill the spaces
    :param align: either left "<", center "^" or right ">"
    :type align: string
    :return: new prepended lines
    :rtype: list
    """
    prepend_pattern = functools.partial("{0:{fill}{align}{indent}}".format, fill=fill, align=align, indent=indent)
    new_lines = []
    for label, line in zip(label_iterable, lines):
        new_lines.append("{}{}".format(prepend_pattern(label), line))
    return new_lines


def indent(string, indent):
    """Indent lines"""
    lines = string.split("\n")
    new_lines = prepend_lines(lines, ['']*len(lines), indent)
    return '\n'.join(new_lines)

# def set_indent(lines, indent):
#     """Reset the indent of lines"""
#     return indent([l.lstrip() for l in lines], indent)
#
#
# def enumerate_lines(lines, indent):
#     """Enumerate lines"""
#     labels = range(len(lines))
#     return prepend_lines(lines, labels, indent)
#
#
# def accumulate_length_of_lines(lines, indent):
#     labels = itertools.accumulate([len(l.strip('\n')) for l in lines], operator.add)
#     return prepend_lines(lines, labels, indent)
#
#
# def accumulate_length_of_first_line(lines, indent):
#     labels = itertools.accumulate([len(l.split('\n')[0].strip('\n')) for l in lines], operator.add)
#     return prepend_lines(lines, labels, indent)


class ViewerAnnotationFlag(object):
    """Flags for annotation directions"""
    FORWARD = ">"
    REVERSE = "<"
    BOTH = "-"


class SequenceRow(object):
    """A row in a :class:`SequenceViewer` instance. Can be comprised of multiple sequences (i.e. lines)
    and can be annotated with 'features'."""

    def __init__(self, lines, labels, indent, start, end):
        """
        SequenceRow constructor

        :param lines: list of lines to display. Lengths of all lines must all be equivalent.
        :type lines: list
        :param labels: list of labels to apply to each line
        :type labels: list
        :param indent: indent to apply to the lines
        :type indent: string
        :param start: start bp of this row
        :type start: int
        :param end: end bp of this row
        :type end: int
        """
        lengths = set([len(r) for r in lines])
        if len(lengths) > 1:
            raise Exception("Cannot format rows that have different lengths")
        self._lines = lines
        self.labels = labels
        self.indent = indent
        self.start = start
        self.end = end
        self.annotations = []

    @property
    def lines(self):
        return prepend_lines(self._lines, self.labels, self.indent)

    @property
    def annotation_lines(self):
        condensed = StringColumn.condense(self.annotations)
        return [str(a.indent(self.indent)) for a in condensed]

    @staticmethod
    def make_annotation(label, span, fill='*'):
        """
        Make an annotation with 'label' spanning inclusive base pairs indices 'span'

        :param label: annotation label
        :type label: basestring
        :param span: the start and end (inclusive) of the annotation
        :type span: tuple
        :param fill: what to fill whitespace with
        :type fill: basestring
        :return:
        :rtype:
        """
        s = ''
        if len(fill) != 1:
            raise Exception("Fill '{}' must be a single character long, not {} characters".format(fill, len(fill)))
        if fill.strip() == '':
            raise Exception("Fill cannot be whitespace")
        sc = StringColumn()
        if len(label) + 1 > span:
            sc.append_string("|<{0:{fill}{align}{indent}}".format(label, fill=' ', align='^', indent=span))
            label = fill * span
        sc.append_string("{0:{fill}{align}{indent}}".format(label, fill=fill, align='^', indent=span))
        return sc

    def absolute_annotate(self, start, end, fill, label):
        """
        Applyt annotation to this row using absolute start and ends for
        THIS row.

        :param start: inclusive start
        :type start: int
        :param end: inclusive end
        :type end: int
        :param fill: what to fill whitespace with
        :type fill: basestring
        :param label: annotation label
        :type label: basestring
        :return: None
        :rtype: None
        """
        span = end - start + 1
        annotation = self.make_annotation(label, span, fill).indent(start)
        self.annotations.append(annotation)

    def annotate(self, start, end, fill, label=''):
        """
        Annotate the sequence row. If 'start' or 'end' is beyond,
        the expected start or end for this row, the annotation will
        automatically be truncated.

        :param start: inclusive start
        :type start: int
        :param end: inclusive end
        :type end: int
        :param fill: what to fill whitespace with
        :type fill:
        :param label: optional label to apply to the annotation
        :type label: basestring
        :return:
        :rtype:
        """
        s = max(start - self.start, 0)
        e = min(end - self.start, len(self)-1)
        return self.absolute_annotate(s, e, fill, label)

    def in_bounds(self, x):
        """
        Checks if the index 'x' is in between row start and end (inclusive)

        :param x: index
        :type x: int
        :return: if in bounds
        :rtype: bool
        """
        return x >= self.start and x <= self.end

    def __len__(self):
        return len(self._lines[0])

    def __str__(self):
        return '\n'.join(self.annotation_lines + self.lines)

#
# class SequenceLabel(object):
#
#     def __init__(self, indent, label=None, pattern=None, indexer=None):
#         self.indent = indent
#         self.index = 0
#         self.label = label
#         if pattern is None:
#             pattern = "{label} {index}"
#         self.pattern = pattern
#         self.indexer = indexer
#
#     def indexers(self):
#         return {
#             "line_length": lambda x: self.index + len(x),
#             "enumerate": lambda x: x + 1
#         }
#
#     def enumerate(self, line):
#         if self.indexer:
#             self.index += self.indexer(line)
#
#     def __str__(self):
#         label = self.patter.format(index=self.index, label=self.label)
#         return "{0:{fill}{align}{indent}".format(label, fill=' ', align='<', indent=self.indent)


class SequenceViewer(object):
    """A class that views longs sets of sequences."""

    class DEFAULTS:
        METADATA_INDENT = 2
        INDENT = 10
        SPACER = '\n'
        WIDTH = 85
        NAME = 'Unnamed'
        DESCRIPTION = ''

    def __init__(self, sequences,
                 sequence_labels=None,
                 indent=DEFAULTS.INDENT,
                 width=DEFAULTS.WIDTH,
                 spacer=DEFAULTS.SPACER,
                 name=DEFAULTS.NAME,
                 description='',
                 metadata=None):
        """
        SequenceViewer constructor

        :param sequences: list of sequences to view
        :type sequences: list
        :param sequence_labels: optional labels to apply to sequence. Include the '{index}' to enumerate the base pairs.
        :type sequence_labels: list
        :param indent: spacing before start of string and start of base pairs
        :type indent: int
        :param width: width of the view window for the sequences (e.g. width=100 would mean rows of at most len 100
                        characters
        :type width: string
        :param spacer: string to apply inbetween rows (default is newline)
        :type spacer: string
        :param name: optional name for this viewer, to be displayed in the header
        :type name: basestring
        :param description: optional description for this viewer
        :type description: basestring
        :param metadata: optional metadata to display in the header
        :type metadata: dict
        """
        assert isinstance(sequences, list)
        seq_lens = set([len(s) for s in sequences])
        if len(seq_lens) > 1:
            raise Exception("Sequence must be same length but found lengths {}".format([len(s) for s in sequences]))

        self._sequences = [str(s) for s in sequences]
        if sequence_labels is None:
            sequence_labels = ["{index}"] + ['']*(len(sequences)-1)
        self._sequence_labels = sequence_labels
        self._indent = indent
        self._width = width
        self._spacer = spacer
        self._rows = self.create_rows()
        if name is None:
            name = self.DEFAULTS.NAME
        self.name = name
        self.metadata = OrderedDict()
        if description:
            self.metadata['Description'] = self.DEFAULTS.DESCRIPTION
        if metadata is not None:
            self.metadata.update(metadata)

    def reset(self):
        self._rows = None

    @property
    def header(self):
        """Return the formatted header and metadata"""
        metadata = '\n'.join('{key}: {val}'.format(key=key, val=val) for key, val in self.metadata.items())
        metadata = indent(metadata, self.DEFAULTS.METADATA_INDENT)
        return "> \"{name}\" ({length}bp)\n{metadata}".format(name=self.name, length=len(self), metadata=metadata)

    @property
    def width(self):
        return self._width

    @width.setter
    def width(self, w):
        self.reset()
        self._width = w

    @property
    def spacer(self):
        return self._spacer

    @spacer.setter
    def spacer(self, s):
        self.reset()
        self._spacer = s

    @property
    def indent(self):
        return self._indent

    @indent.setter
    def indent(self, i):
        self.reset()
        self._indent = i

    @property
    def sequences(self):
        return self._sequences

    @property
    def rows(self):
        if not self._rows:
            self._rows = self.create_rows()
        return self._rows

    def create_rows(self):
        """Create :class:`SequenceRow` instances for the set of sequences"""
        lines = []
        for seq in self.sequences:
            line = to_lines(str(seq), width=self.width)
            lines.append(line)
        interleafed = functools.reduce(lambda x, y: x + y, zip(*lines))
        chunks = chunkify(interleafed, len(self.sequences))
        rows = []
        index = 0
        for chunk in chunks:
            labels = [str(l).format(index=index) for l in self._sequence_labels]
            rows.append(SequenceRow(chunk, labels, self.indent, index, min(index + self.width - 1, len(self))))
            index += len(chunk[0])
        return rows

    def annotate(self, start, end, label=None, direction=None):
        """
        Annotates this viewer object starting from 'start' to 'end' inclusively.

        :param start: inclusive start
        :type start: int
        :param end: inclusive end
        :type end: int
        :param label: optional label to apply to the annotation
        :type label: basestring
        :param direction: the direction of the annotation ('<', '>', '^') to fill in whitespace
        :type direction: string
        :return: None
        :rtype: None
        """
        if direction is None:
            direction = ViewerAnnotationFlag.BOTH
        if label is None:
            label = ''
        for row in self.rows:
            if end >= row.start and start <= row.end:
                row.annotate(start, end, label=label, fill=str(direction))

    def print(self):
        print(str(self))

    def __len__(self):
        return len(self.sequences[0])

    def __str__(self):
        spacer = self.spacer
        if spacer is None:
            spacer = ''
        s = "{header}\n\n".format(header=self.header)
        s += '\n{}'.format(spacer).join([str(r) for r in self.rows])
        return s





