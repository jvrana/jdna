import pytest

from jdna.reaction import Reaction
from jdna.sequence import Sequence

@pytest.fixture(scope='function')
def template():
    return Sequence.random(500)

@pytest.mark.parametrize('o1', [0])
@pytest.mark.parametrize('o2', [0])
@pytest.mark.parametrize('l1', [15])
@pytest.mark.parametrize('l2', [15])
def test_pcr(template, l1, l2, o1, o2):
    pos1 = 100
    pos2 = 200
    p1 = Sequence('N'*o1) + template[pos1:pos1+l1]
    p2 = Sequence('N'*o2) + template[pos2-l2:pos2].reverse_complement()

    products = Reaction.pcr(template, p1, p2)
    assert len(products) == 1

