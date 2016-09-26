from collections import defaultdict
from hydradna.nucleotide import Nucleotide
from hydradna.linkedset import LinkedSet

class Sequence(LinkedSet):

    def __init__(self, first_nt=None, sequence=None, name='unknown'):
        super(Sequence, self).__init__(first=first_nt, data_sequence=sequence)
        self.name = name
        self.description = ''

    def initialize(self, data_sequence):
        self.first = Nucleotide(data_sequence[0])
        current = self.first
        for d in data_sequence[1:]:
            new = Nucleotide(d)
            current.set_next(new)
            current = new

    def add_feature(self, i, j, feature, start_index=0):
        s = self.get()
        start = s[i]
        end = s[j]
        feature_nts = start.fwd(stop_link=end)
        if not feature_nts[-1] == end:
            raise IndexError("Feature index error")
        for i, n in enumerate(feature_nts):
            n.add_feature(feature, i+start_index)
        feature.length = len(feature_nts)-1
        return feature

    def _features_to_i(self):
        features = defaultdict(list)
        for i, x in enumerate(self.get()):
            for f in x.features:
                features[f].append((x,i))
        return features

    def get_features(self):
        features_to_nts = self._features_to_i()
        feature_info = {}
        for feature in features_to_nts:
            nt_to_i_list = features_to_nts[feature]
            nts, indices = zip(*nt_to_i_list)
            first = nts[0].feature_rev(feature)[-1]
            last = nts[0].feature_fwd(feature)[-1]
            nt_to_i = dict(zip(self.get(), range(len(self))))
            feature_range = (first.features[feature], last.features[feature])
            pos_ranges = (nt_to_i[first], nt_to_i[last])
            feature_info[feature] = [pos_ranges, feature_range]
        return feature_info

    def print_features(self):
        features = self.get_feature_pos()
        print "Features:",
        for f in features:
            e = f.end
            if e is None:
                e = ''
            span = '[{}...{}]'.format(f.start, e)
            if f.end == None and f.start == 0:
                span = ''
            print '{}{} {}'.format(f.name, span, f.type),
            ranges = list(Utilities.group_ranges(features[f]))
            feature_ranges = []
            for r in ranges:
                print '{}...{},'.format(r[0], r[-1]),
        print

    def find_feature(self, name=None):
        found = []
        for feature in self.get_features():
            if feature.name == name:
                found.append(feature)
        return found

    def create_feature(self, name, type, start, end):
        f = Feature(name, type)
        self.add_feature(start, end, f)
        return f

    def complement(self):
        curr = self.get_first()
        visited = set()
        while curr and curr not in visited:
            visited.add(curr)
            curr.to_complement()
            curr = curr.next()
        return self

    def reverse_complement(self):
        self.reverse()
        self.complement()
        return self

    def _chop(self, i, cut_method):
        self._inbounds(i)
        nt = self.get()[i]
        n = cut_method(nt)
        if n is None:
            return self, None
        remaining = Sequence(first_nt=n)
        return self, remaining

    def chop_prev(self, i):
        self._chop(i, lambda x: x.cut_prev())

    def chop_prev(self, i):
        self._chop(i, lambda x: x.cut_next())

    def cut(self, i, cut_prev=True):
        fragments = super(Sequence, self).cut(i, cut_prev)
        fragments = [Sequence(first_nt=f.get_first()) for f in fragments]
        return fragments

    def fuse(self, seq):
        f = self.get_first().find_last()
        l = seq.get_first()
        f.set_next(l)
        return self

    def __add__(self, other):
        return deepcopy(self).fuse(deepcopy(other))

    def __repr__(self):
        return str(self)