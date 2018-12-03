import pytest

from jdna.reaction import Reaction
from jdna.sequence import Sequence


@pytest.mark.parametrize('o1', [10, 0])
@pytest.mark.parametrize('o2', [10, 0])
@pytest.mark.parametrize('primer1_inset', [
    pytest.param(30, id="inset1=30bp"),
    pytest.param(15, id="inset1=15bp"),
])
@pytest.mark.parametrize('primer2_inset', [
    pytest.param(30, id="inset2=30bp"),
    pytest.param(15, id="inset2=15bp"),
])
@pytest.mark.parametrize('circular_template', [False, True])
def test_pcr(primer1_inset, primer2_inset, o1, o2, circular_template):
    template = Sequence.random(100)

    primer_len = 20
    p1_start = primer1_inset
    p1_end = primer1_inset + primer_len

    p2_end = len(template) - primer2_inset
    p2_start = p2_end - primer_len

    print(p1_start)
    print(p1_end)
    print(p2_start)
    print(p2_end)

    p1 = Sequence('N'*o1) + template[p1_start:p1_end]
    p2 = Sequence('N'*o2) + template[p2_start:p2_end].reverse_complement()

    template.circularize()
    if circular_template:
        template.reindex(int((len(template) - primer2_inset - primer1_inset)/2))

    products = Reaction.pcr(template, [p1, p2])
    assert len(products) == 1
    expected_len = len(template) - primer1_inset - primer2_inset + o1 + o2
    assert len(products[0]) == expected_len, "Product should be {} long".format(expected_len)


def test_pcr_expected_sequence():
    anneal1 = Sequence.random(20)
    anneal2 = Sequence.random(20)
    overhang1 = Sequence.random(20)
    overhang2 = Sequence.random(20)
    buffer = Sequence('N'*100)

    template = buffer + anneal1 + buffer + anneal2 + buffer
    primer1 = overhang1 + anneal1
    primer2 = overhang2 + anneal2.copy().rc()

    products = Reaction.pcr(template, [primer1, primer2])
    assert len(products) == 1

    print(anneal2)
    print(anneal2.copy().reverse_complement())
    print(overhang2)
    print(overhang2.copy().reverse_complement())

    assert str(products[0]) == str(overhang1 + anneal1 + buffer + anneal2 + overhang2.copy().reverse_complement())