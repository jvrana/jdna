from hydradna.models.linkedset import *
from hydradna.models.feature import *


class Nucleotide(Link):

    base_pairing = dict(zip(
        ['a', 't', 'c', 'g', 'A', 'T', 'C', 'G'],
        ['t', 'a', 'g', 'c', 'T', 'A', 'G', 'C']
    ))

    def __init__(self, base):
        super(Nucleotide, self).__init__(base)
        self.features = set()
        self.pair = None

    def base(self):
        return self.data

    def to_complement(self):
        self.data = self.base_pairing[self.data]

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
            # If this feature doesn't span, then continue
            if f not in x2.features:
                continue

            # Grab the sequences for the split feature
            frag1 = x1._feature_rev(f)
            frag2 = x2._feature_fwd(f)

            # Make two copies of the feature
            f1 = deepcopy(f)
            f2 = deepcopy(f)

            # modify end of f1
            if f1.end is None:
                f1.end = 0
            f1.end += f1.start + len(frag1) - 1

            # modify start of f2
            if f2.start is None:
                f2.start = 0
            f2.start += len(frag1)
            # leave end of f2 alone

            try:
                for n in frag1:
                    n.features.remove(f)
                    n.features.add(f1)
            except:
                pass
            try:
                for n in frag2:
                    n.features.remove(f)
                    n.features.add(f2)
            except:
                pass

    def _cut(self, cut_prev=True):
        for f in self.features:
            self.split_features(split_prev=cut_prev)
        nxt = None
        if cut_prev:
            nxt = super(Nucleotide, self).cut_prev()
        else:
            nxt = super(Nucleotide, self).cut_next()
        return nxt

    def cut_prev(self):
        return self._cut(cut_prev=True)

    def cut_next(self):
        return self._cut(cut_prev=False)

    def _feature_fwd(self, feature):
        stop = lambda x: feature not in x.features
        return self._propogate(lambda x: x.next(), stop_criteria=stop)

    def _feature_rev(self, feature):
        stop = lambda x: feature not in x.features
        return self._propogate(lambda x: x.prev(), stop_criteria=stop)

    def replace_feature(self, old_feature, new_feature):
        self.features.remove(old_feature)
        self.features.add(new_feature)

    @staticmethod
    def fuse_features(n1, n2):
        del1, del2, fused_features = set(), set(), set()
        if n1 is None:
            return
        if n2 is None:
            return
        for f1 in n1.features:
            for f2 in n2.features:
                if f1.name == f2.name:
                    if not f1.end is None and f1.end + 1 == f2.start:

                        # Setup features to delete
                        del1.add(f1)
                        del2.add(f2)

                        # Create fused feature
                        f = Feature(f1.name, f1.type)
                        f.start = f1.start
                        f.end = f2.end

                        # Setup features
                        fused_features.add(f)

        for d in del1:
            for n in n1._feature_rev(d):
                n.features.remove(d)
                [n.features.add(x) for x in fused_features]
        for d in del2:
            for n in n2._feature_fwd(d):
                n.features.remove(d)
                [n.features.add(x) for x in fused_features]

    def set_next(self, nucleotide):
        super(Nucleotide, self).set_next(nucleotide)
        Nucleotide.fuse_features(self, nucleotide)

    def set_prev(self, nucleotide):
        super(Nucleotide, self).set_prev(nucleotide)
        Nucleotide.fuse_features(nucleotide, self)

    def equivalent(self, other):
        return str(self.data).upper() == str(other.data).upper()