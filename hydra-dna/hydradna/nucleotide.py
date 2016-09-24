from hydradna.link import Link
from itertools import combinations

class Nucleotide(Link):

    base_pairing = dict(zip(
        ['a', 't', 'c', 'g', 'A', 'T', 'C', 'G'],
        ['t', 'a', 'g', 'c', 'T', 'A', 'G', 'C']
    ))

    def __init__(self, base):
        super(Nucleotide, self).__init__(base)
        self.features = {}

    def base(self):
        return self.data

    def equivalent(self, other):
        return str(self.base()).upper() == str(other.base()).upper()

    def to_complement(self):
        self.data = self.base_pairing[self.data]

    def set_next(self, nucleotide):
        super(Nucleotide, self).set_next(nucleotide)
        Nucleotide.fuse_features(self, nucleotide)

    def set_prev(self, nucleotide):
        super(Nucleotide, self).set_prev(nucleotide)
        Nucleotide.fuse_features(nucleotide, self)

    def cut_prev(self):
        return self._cut(cut_prev=True)

    def cut_next(self):
        return self._cut(cut_prev=False)

    def _cut(self, cut_prev=True):
        for f in self.features:
            self.split_features(split_prev=cut_prev)
        nxt = None
        if cut_prev:
            nxt = super(Nucleotide, self).cut_prev()
        else:
            nxt = super(Nucleotide, self).cut_next()
        return nxt

    def _feature_fwd(self, feature):
        stop = lambda x: feature not in x.features
        return self._propogate(lambda x: x.next(), stop_criteria=stop)

    def _feature_rev(self, feature):
        stop = lambda x: feature not in x.features
        return self._propogate(lambda x: x.prev(), stop_criteria=stop)

    def add_feature(self, feature, pos):
        self.features[feature] = pos
        return feature

    def replace_feature(self, old_feature, new_feature):
        self.features[new_feature] = self.features[old_feature]
        self.remove_feature(old_feature)

    def copy_features_from(self, other):
        for f in other.features:
            i = other.features[f]
            if f not in self.features:
                self.add_feature(f, i)
        self._remove_overlapping_features()

    def remove_feature(self, feature):
        del self.features[feature]

    def get_feature_span(self, feature):
        start = self._feature_rev(feature)[-1]
        end = self._feature_fwd(feature)[-1]
        return (start.features[feature], end.features[feature])

    # def update_feature_span(self, feature, delta_i):
    #     start = self._feature_rev(feature)[-1]
    #     for n in start._feature_fwd(feature):
    #         n.features[feature] += delta_i

    def _remove_overlapping_features(self):
        # type: () -> Nucleotide
        feature_pairs = combinations(self.features.keys(), 2)
        tobedel = set()
        for f1, f2 in feature_pairs:
            if f1.name == f2.name:
                tobedel.add(f2)
        for tob in tobedel:
            self.remove_feature(tob)

    @staticmethod
    def fuse_features(n1, n2):
        if n1 is None:
            return
        if n2 is None:
            return

        delset = set()

        for f1 in n1.features:
            for f2 in n2.features:
                f1_pos = n1.features[f1]
                f2_pos = n2.features[f2]
                f1_copy = deepcopy(f1)
                # same name & consecutive position
                if f1 is f2:
                    continue
                if f1.name == f2.name and f1_pos + 1 == f2_pos:
                    delset.add((f1, f2, f1_copy))
        for f1, f2, f1_copy in delset:
            for n in n1._feature_rev(f1):
                try:
                    n.replace_feature(f1, f1_copy)
                except:
                    pass
            for n in n2._feature_fwd(f2):
                try:
                    n.replace_feature(f2, f1_copy)
                except:
                    pass


    def split_features(self, split_prev=True):
        x1 = self.prev()
        x2 = self
        if not split_prev:
            # then split_next
            x1 = self
            x2 = self.next()
        # If at the end, no splitting is necessary
        if x1 is None or x2 is None:
            return
        for f in x1.features:
            # If this feature spans
            if f in x2.features:
                # Grab the sequences for the split feature
                frag1 = x1._feature_rev(f)
                frag2 = x2._feature_fwd(f)

                # Make two copies of the feature
                f1 = deepcopy(f)
                f2 = deepcopy(f)

                # Swap original feature for copy
                for n in frag1:
                    n.replace_feature(f, f1)
                for n in frag2:
                    n.replace_feature(f, f2)