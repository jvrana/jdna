'''
Project: jdna
File: test_plasmid_viewer
Author: Justin
Date: 2/13/17

Description: 

'''

import coral
import os
import json
seq = coral.seqio.read_dna("designs/pmodkan-ho-pact1-z4-er-vpr.gb")
print seq.features[0].__dict__
print seq.__dict__.keys()

def to_json(c):
    d = c.__dict__
    d['width'] = 700
    d['height'] = d['width']
    d['length'] = len(c)
    d['features'] = [x.__dict__ for x in c.features]
    del d['top']
    del d['bottom']
    print json.dumps(d)

to_json(seq)