from hydradna.models.linkedset import *
from hydradna.models.nucleotide import *
from hydradna.utilities import Utilities
from collections import defaultdict


class Sequence(LinkedSet):

    def __init__(self, first_nt=None, sequence=None, name=None):
        super(Sequence, self).__init__(first=first_nt, data_sequence=sequence)
        self.name = name

    def initialize(self, data_sequence):
        self.first = Nucleotide(data_sequence[0])
        current = self.first
        for d in data_sequence[1:]:
            new = Nucleotide(d)
            current.set_next(new)
            current = new

    def add_feature(self, i, j, feature):
        s = self.get()
        start = s[i]
        end = s[j]
        feature_nts = start.fwd(stop_link=end)
        if not feature_nts[-1] == end:
            raise IndexError("Feature index error")
        for n in feature_nts:
            n.features.add(feature)

    def get_features(self):
        features = defaultdict(list)
        for i, x in enumerate(self.get()):
            for f in x.features:
                features[f].append(i)
        return features

    def get_feature_ranges(self):
        features = self.get_features()
        for f in features:
            ranges = list(Utilities.group_ranges(features[f]))
            feature_ranges = []
            for r in ranges:
                feature_ranges.append([r[0], r[-1]])
            features[f] = feature_ranges
        return features

    def print_features(self):
        features = self.get_feature_ranges()
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
        found = set()
        for feature in self.get_features():
            if feature.name == name:
                found.add(feature)
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

    def cut(self, i, cut_prev=True):
        fragments = super(Sequence, self).cut(i, cut_prev)
        fragments = [Sequence(first_nt=f.get_first()) for f in fragments]
        return fragments

    def fuse(self, seq):
        f = self.get_first().find_last()
        l = seq.get_first()
        f.set_next(l)
        return self

    def __repr__(self):
        return str(self)