import functools
import itertools
import operator


class NumberStrategy:
    NUM = "line_number"
    LENGTH = "line_length"


def chunks(iterable, n):
    chunk = None
    for i, x in enumerate(iterable):
        if i % n == 0:
            if chunk is not None:
                yield chunk
            chunk = []
        chunk.append(x)
    yield chunk


def to_lines(string, width):
    """Converts a string to lines of length <= width."""
    lines = []
    for i in range(0, len(string), width):
        lines.append(string[i : i + width])
    return lines


def rows_to_chunks(rows, width=75):
    if not rows:
        return []
    lengths = {len(r) for r in rows}
    if len(lengths) > 1:
        raise Exception("Cannot format rows that have different lengths")

    lines = []
    for row in rows:
        line = to_lines(row, width=width)
        lines.append(line)
    #     print(lines)
    interleafed = functools.reduce(lambda x, y: x + y, zip(*lines))
    return list(chunks(interleafed, len(rows)))


def prepend_lines(lines, label_iterable, indent, fill=" ", align="<"):
    """

    :param lines:
    :type lines:
    :param indent:
    :type indent: int
    :param fill: default ' '
    :type fill: what to fill the spaces
    :param align: either left "<", center "^" or right ">"
    :type align: string
    :return:
    :rtype:
    """
    prepend_pattern = functools.partial(
        "{0:{fill}{align}{indent}}".format, fill=fill, align=align, indent=indent
    )
    new_lines = []
    for label, line in zip(label_iterable, lines):
        print(line)
        new_lines.append("{}{}".format(prepend_pattern(label), line))
    return new_lines


def indent(lines, indent):
    return prepend_lines(lines, [""] * len(lines), indent)


def set_indent(lines, indent):
    return indent([l.lstrip() for l in lines], indent)


def enumerate_lines(lines, indent):
    labels = range(len(lines))
    return prepend_lines(lines, labels, indent)


def accumulate_length_of_lines(lines, indent):
    labels = itertools.accumulate([len(l.strip("\n")) for l in lines], operator.add)
    return prepend_lines(lines, labels, indent)


def accumulate_length_of_first_line(lines, indent):
    labels = itertools.accumulate(
        [len(l.split("\n")[0].strip("\n")) for l in lines], operator.add
    )
    return prepend_lines(lines, labels, indent)


def number_lines(lines, indent=None, numbering_strategy=NumberStrategy.LENGTH, step=1):
    if indent is None:
        indent = 10
    index = 0
    indices = []
    for i, line in enumerate(lines):
        if i % step == 0:
            indices.append(index)
            if numbering_strategy == NumberStrategy.LENGTH:
                index += len(line)
            elif numbering_strategy == NumberStrategy.NUM:
                index += 1
        else:
            indices.append("")

    new_lines = []
    for index, line in zip(indices, lines):
        index_str = "{0:{fill}{align}{indent}}".format(
            index, fill=" ", align="<", indent=indent
        )
        new_lines.append("{}{}".format(index_str, line))
    return new_lines


def interleaf(rows, width=75, spacer=None, number=False, indent=None):
    """Interleaf lines that have the same lengths."""
    if not rows:
        return ""
    lengths = {len(r) for r in rows}
    if len(lengths) > 1:
        raise Exception("Cannot format rows that have different lengths")

    lines = []
    for row in rows:
        line = to_lines(row, width=width)
        lines.append(line)
    interleafed = functools.reduce(lambda x, y: x + y, zip(*lines))

    if number:
        interleafed = number_lines(interleafed, step=len(rows), indent=indent)

    with_spacer = []
    for i, l in enumerate(interleafed):
        if spacer is not None and i % len(rows) == 0:
            with_spacer.append(spacer)
        with_spacer.append(l)

    return "\n".join(with_spacer)


def group_by_line_lengths(lines):
    # group by line len
    groups = [[]]
    length = len(lines[0])
    for line in lines:
        if len(line) != length:
            groups.append([])
        if line:
            groups[-1].append(line)
    return groups


def format_sequence(s, width=75, spacer=""):
    groups = group_by_line_lengths(s.split("\n"))
    return "\n".join(
        [interleaf(g, width=width, spacer=spacer, number=True) for g in groups]
    )
