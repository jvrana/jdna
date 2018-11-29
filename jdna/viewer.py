import functools
import itertools
import operator
from enum import Enum


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


class AnnotationFlag(object):
    FORWARD = ">"
    REVERSE = "<"
    BOTH = "^"


class SequenceRow(object):

    def __init__(self, lines, labels, indent, start, end):
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
        annotations = []
        for a in self.annotations:
            annotations.append(indent(a, self.indent))
        return annotations

    @staticmethod
    def make_annotation(label, span, fill='*'):
        s = ''
        if len(label) > span:
            s += "|<{0:{fill}{align}{indent}}\n".format(label, fill=' ', align='^', indent=span)
            label = fill * span
        s += "{0:{fill}{align}{indent}}".format(label, fill=fill, align='^', indent=span)
        return s

    def absolute_annotate(self, start, end, fill, label):
        span = end - start + 1
        annotation = self.make_annotation(label, span, fill)
        annotation_lines = [' '*start + a for a in annotation.split('\n')]
        self.annotations.append(
            '\n'.join(annotation_lines)
        )

    def annotate(self, start, end, fill, label=''):
        s = max(start - self.start, 0)
        e = min(end - self.start, len(self)-1)
        return self.absolute_annotate(s, e, fill, label)

    def in_bounds(self, x):
        return x >= self.start and x <= self.end

    def __len__(self):
        return len(self._lines[0])

    def __str__(self):
        return '\n'.join(self.annotation_lines + self.lines)


class SequenceViewer(object):

    def __init__(self, sequences, indent=10, width=85, spacer='\n', name=None):
        assert isinstance(sequences, list)
        seq_lens = set([len(s) for s in sequences])
        if len(seq_lens) > 1:
            raise Exception("Sequence must be same length but found lengths {}".format([len(s) for s in sequences]))
        self._sequences = [str(s) for s in sequences]
        self._indent = indent
        self._width = width
        self._spacer = spacer
        self._rows = self.create_rows()
        if name is None:
            name = 'Unnamed'
        self.name = name

    def reset(self):
        self._rows = None

    @property
    def header(self):
        return "> \"{name}\" ({length}bp)".format(name=self.name, length=len(self))

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
        lines = []
        for seq in self.sequences:
            line = to_lines(str(seq), width=self.width)
            lines.append(line)
        interleafed = functools.reduce(lambda x, y: x + y, zip(*lines))
        chunks = chunkify(interleafed, len(self.sequences))
        rows = []
        index = 0
        for chunk in chunks:
            labels = [index] + [''] * (len(self.sequences) - 1)
            rows.append(SequenceRow(chunk, labels, self.indent, index, min(index + self.width - 1, len(self))))
            index += len(chunk[0])
        return rows

    def annotate(self, start, end, label=None, direction=None):
        if direction is None:
            direction = AnnotationFlag.BOTH
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





