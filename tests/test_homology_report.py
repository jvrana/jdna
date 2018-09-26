'''
Project: jdna
File: test_homology_report.py
Author: Justin
Date: 1/12/17

Description: 

'''

import glob
import json
from jdna.core import Sequence, Reaction, Feature, Convert
from copy import copy


def test_long_gibson_assembly():
    dnas = []
    for file in glob.glob('test_data/Frag*json'):
        with open(file, 'r') as f:
            bseq = json.load(f)
            dnas.append(Convert.from_benchling(bseq))
    dnas_copy = [copy(x) for x in dnas]
    h_report = Reaction.homology_report(dnas)
    Reaction.print_homology_report(h_report)
    products = Reaction.homology_assembly(dnas, True)
    Reaction.homology_report(dnas)

    # Verify original dnas did not change
    assert len(products) == 1
    assert len(products[0]) > 1000

    assert [len(x) for x in dnas] == [len(x) for x in dnas_copy]
    for i in range(len(dnas)):
        assert str(dnas[i]) == str(dnas_copy[i])

    assert len(dnas) == len(dnas_copy)