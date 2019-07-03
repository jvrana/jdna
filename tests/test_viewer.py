import pytest
from jdna import Sequence
from jdna.viewer import SequenceViewer, ViewerAnnotationFlag


def test_viewer():
    s = Sequence.random(500)
    viewer = SequenceViewer([s, s.copy().complement()])
    print(viewer.rows)
    viewer.print()


@pytest.mark.parametrize(
    "length,width,expected_num_rows",
    [(1000, 500, 2), (500, 1000, 1), (1000, 100, 10), (1001, 100, 11)],
)
def test_viewer_expected_rows(length, width, expected_num_rows):
    viewer = SequenceViewer([Sequence.random(length)], width=width)
    assert len(viewer.rows) == expected_num_rows


def test_viewer_print():
    SequenceViewer([Sequence.random(500)]).print()


@pytest.mark.parametrize(
    "positions,expected_num_annotations",
    [
        ([(0, 50), (75, 99), (150, 199)], 3),
        ([(0, 50)], 1),
        ([(75, 99), (75, 99)], 2),
        ([(75, 99), (75, 100)], 3),
    ],
)
def test_num_annotations(positions, expected_num_annotations):
    viewer = SequenceViewer([Sequence.random(1000)], width=100)
    for p in positions:
        viewer.annotate(*p)

    all_annotations = []
    for row in viewer.rows:
        all_annotations += row.annotations
    assert len(all_annotations) == expected_num_annotations
    viewer.print()


@pytest.mark.parametrize(
    "start,end,rows,direction,fill",
    [
        (50, 70, [0], ViewerAnnotationFlag.BOTH, "-"),
        (50, 70, [0], None, "-"),
        (80, 99, [0], ViewerAnnotationFlag.BOTH, "-"),
        pytest.param(
            80,
            100,
            [0, 1],
            ViewerAnnotationFlag.BOTH,
            "-",
            id="annotation spans 2 rows",
        ),
        pytest.param(
            80,
            200,
            [0, 1, 2],
            ViewerAnnotationFlag.BOTH,
            "-",
            id="annotation spans 3 rows",
        ),
        (50, 70, [0], ViewerAnnotationFlag.FORWARD, ">"),
        (50, 70, [0], ViewerAnnotationFlag.REVERSE, "<"),
        (500, 699, [], ViewerAnnotationFlag.REVERSE, "<"),
    ],
)
def test_annotate_viewer_with_fill(start, end, rows, direction, fill):
    viewer = SequenceViewer([Sequence.random(500)], width=100)
    viewer.annotate(start, end, fill=direction)
    viewer.name = "{}->{}".format(start, end)
    viewer.print()
    if rows:
        assert fill in str(viewer)
    rows_with_fill = []
    for i, r in enumerate(rows):
        if fill in str(viewer.rows[r]):
            rows_with_fill.append(r)
    assert rows_with_fill == rows
    viewer.print()


@pytest.mark.parametrize(
    "start,end,label,expected_labels,no_labels",
    [
        (0, 50, "mylabel", ["mylabel", "-"], []),
        (0, 2, "mylabel", ["-", "mylabel"], []),
        (0, 5, "mylabel", ["mylabel"], []),
        (0, 7, "mylabel", ["mylabel"], []),
        (2, 6, "mylabel", ["mylabel"], []),
    ],
)
def test_annotate_viewer_with_label(start, end, label, expected_labels, no_labels):
    viewer = SequenceViewer([Sequence.random(500)], width=100)
    viewer.annotate(start, end, label=label, background="blue")
    for expected_label in expected_labels:
        assert expected_label in str(viewer.rows[0])
    for no_label in no_labels:
        assert no_label not in str(viewer.rows[0])
    viewer.print()
