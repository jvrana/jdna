

from jdna.core import Sequence


seq = Sequence(sequence='agttggagcg')
seq.make_cyclic()
seq.create_feature('feat', 't', 5, 4)
print seq.get()[5].features
print seq.get_features()
print seq.get()[5].features
# print seq.get_first().features
# print seq.get_features()
# seq.linearize()
# print seq.get_features()