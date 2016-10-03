from jdna.core import Reaction, Convert
from copy import copy
import json
import glob


def test_long_gibson_assembly():
    dnas = []
    for file in glob.glob('test_data/Fragment*json'):
        with open(file, 'r') as f:
            bseq = json.load(f)
            dnas.append(Convert.from_benchling(bseq))
    dnas_copy = [copy(x) for x in dnas]
    Reaction.evaluate_assembly(dnas)
    products = Reaction.cyclic_assembly(dnas)
    Reaction.evaluate_assembly(dnas)
    # Verify original dnas did not change
    assert len(products) == 1
    assert len(products[0]) > 1000

    assert [len(x) for x in dnas] == [len(x) for x in dnas_copy]
    for i in range(len(dnas)):
        assert str(dnas[i]) == str(dnas_copy[i])

    assert len(dnas) == len(dnas_copy)
