import pytest
from jdna.format import interleaf, format_sequence
from jdna import format, Sequence
from jdna.viewer import SequenceViewer


@pytest.mark.parametrize('spacer', [None, ''])
@pytest.mark.parametrize('width', [50, 88])
def test_interleaf(spacer, width):
    a = 'a'*960
    b = 'b'*960
    s = interleaf('\n'.join([a, b]), spacer=spacer, width=width, number=True)
    print(s)


@pytest.mark.parametrize('spacer', [None, ''])
@pytest.mark.parametrize('width', [50, 88])
def test_format(spacer, width):
    a = 'a'*960
    b = 'b'*960
    c = 'c'*10
    s = format_sequence('\n'.join([a, b, c]), spacer=spacer, width=width)
    print(s)


def test_prepend():
    s = 'prepend'
    lines = [s, s, s]
    print(format.accumulate_length_of_lines(lines, 10))


def test_viewer():
    s = Sequence.random(500)
    viewer = SequenceViewer([s, s.copy().complement()])
    print(viewer.rows)
    viewer.print()