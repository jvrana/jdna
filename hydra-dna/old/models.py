import collections
from itertools import groupby
from copy import deepcopy
import warnings
import uuid

class HydraSequenceException(Exception):
    """ Generic exception for HydraSequence errors """


class Nucleotide(object):
    """
    A nucleotide represents a single basepair in a DNA or RNA sequence.
    Each nucleotide can be connected to another nucleotide through the _next
    and _prev pointers (representing backbone bonds). Nuceotide pairs are connected
    through the _pair pointer. Each nucelotide
    """

    def __init__(self, base):
        self.base = base
        self.features = set()
        self.__next = None
        self.__prev = None
        self.__pair = None

    def get_next(self):
        return self.__next

    def get_prev(self):
        return self.__prev

    def add_next(self, base_letter):
        new_nt = Nucleotide(base_letter)
        self._set_next(new_nt)
        return new_nt

    def cut_next(self):
        nt = self.get_next()
        if nt is not None:
            nt.__assign_prev(None)
        self.__assign_next(None)
        return nt

    def cut_prev(self):
        pr = self.get_prev()
        if pr is not None:
            pr.__assign_next(None)
        self.__assign_prev(None)
        return pr

    def remove_next(self):
        nt = self.get_next()
        self._set_next(nt.next())
        del nt

    def remove_prev(self):
        nt = self.get_prev()
        self._set_prev(nt.prev())
        del nt

    #TODO: Fix pair methods
    def _set_pair(self, nt):
        nt.__pair = self
        self.__pair = nt

    def __assign_next(self, nt):
        '''
        This is the only method that should ever touch __next
        :param nt:
        :return:
        '''
        self.__next = nt

    def __assign_prev(self, nt):
        '''
        This is the only method that should ever touch __prev
        :param nt:
        :return:
        '''
        self.__prev = nt

    def _set_next(self, nt):
        '''
        This is the only method that should ever call __set
        :param nt:
        :return:
        '''
        nt.__assign_prev(self)
        self.__assign_next(nt)

    def _set_prev(self, nt):
        nt.set_next(self)

    def _propogate(self, next_method, stop=None, stop_criteria=None):
        visited = []
        next = self
        while next is not None:
            if stop_criteria is not None and stop_criteria(next):
                break
            visited.append(next)
            if next is stop:
                break
            next = next_method(next)
            if next is visited[0]:
                break
        return visited

    def fwd(self, stop_nt=None, stop_criteria=None):
        return self._propogate(
            lambda x: x.next(),
            stop=stop_nt,
            stop_criteria=stop_criteria)

    def rev(self, stop_nt=None, stop_criteria=None):
        return self._propogate(
            lambda x: x.prev(),
            stop=stop_nt,
            stop_criteria=stop_criteria)

    def find_first(self):
        seq = self.rev()
        return seq[-1]

    def find_last(self):
        seq = self.fwd()
        return seq[-1]

    def __repr__(self):
        return self.base

class Feature(object):
    def __init__(self, name, feature_type):
        self.name = name
        self.type = feature_type
        self.splice_start = [] #records the indicies of the original feature that exists after splicing
        self.splice_end = []
        self.uuid = uuid.uuid4()

    def __repr__(self):
        return "<FeatureName: {}, Type: {}".format(self.name, self.type)


# class dsdna(Sequence):
#     pass

# class rna(Sequence):


class Sequence(object):
    """
    A sequence is a series of linked base pairs. You access the sequence only
    from the first base pair and recursively accessing the next base pair. In
    this way, each base pair can maintain its topologies through manipulations.

    Each sequence MUST have unique nucleotide constructors for each basepair unless
    you explicitly wish to manipulate base pairs across different sequences (e.g. chaning a common part...)

    self.first_base is always the index=0
    """

    def __init__(self, name='', first_nucleotide=None, sequence=None, circular=False):
        self.name = name
        if first_nucleotide is not None:
            assert isinstance(first_nucleotide, Nucleotide)
            self.first_base = first_nucleotide
            # self.reindex()
        if sequence is not None:
            assert isinstance(sequence, basestring)
            self.initialize(sequence)
        if circular:
            self.circularize()

    def circularize(self):
        if self.is_circular():
            warnings.warn("Sequence already circular!")
            return
        firstnt = self.first_base.find_first()
        lastnt = self.first_base.find_last()
        firstnt.set_prev(lastnt)

    def reindex(self, index):
        self._inbounds(index)
        if not self.is_circular():
            raise HydraSequenceException("Cannot reindex a linear sequence")
        seq = self.seq()
        self.first_base = seq[index]

    def linearize(self, index=0):
        if self.is_circular():
            self.first_base = self.seq()[index]
            self.first_base.cut_prev()
        else:
            warnings.warn('Sequence is already linear')

    def is_circular(self):
        last = self.first_base.find_last()
        if last.next() == self.first_base and self.first_base.get_prev() == last:
            return True
        return False

    def initialize(self, sequence):
        current_base = Nucleotide(sequence[0])
        self.first_base = current_base
        self._sequence_text = sequence
        for n in sequence[1:]:
            current_base = current_base.add_next(n)

    def seq(self, start_nt=None, end_nt=None):
        if start_nt is None:
            start_nt = self.first_base
        return start_nt.fwd(stop_nt=end_nt)

    def add_feature(self, feature, start, end, overwrite=False):
        existing_features = self.features()
        if feature in existing_features:
            if overwrite:
                self.remove_feature(feature)
            else:
                raise HydraSequenceException("Feature {} exists".format(feature))
        sequence = self.seq()

        self._inbounds_se(start, end)
        for nt in self.get(start,end):
            if feature not in nt.features:
                nt.features.add(feature)

    def remove_feature(self, feature):
        sequence = self.seq()
        for nt in sequence:
            if feature in nt.features:
                nt.features.remove(feature)

    def nt_to_i(self, nt):
        seq = self.seq()
        indices = range(len(seq))
        return dict(zip(seq, indices))[nt]

    def _get_feature_sequence(self, feature, nt):
        feature_sequence = nt.rev(stop_criteria=lambda x: feature not in x.features)[::-1][:-1] + nt.fwd(
            stop_criteria=lambda x: feature not in x.features)

    def features(self):
        """
        Gets a dictionary of features to locations

        """
        seq = self.seq()
        features = collections.defaultdict(list)
        for index, nt in enumerate(seq):
            for f in nt.features:
                features[f].append(index)
        for f in features:
            ranges = list(self._ranges(features[f]))
            if len(ranges) == 1:
                features[f] = (ranges[0][0], ranges[0][-1])
            elif len(ranges) == 2:
                if self.is_circular():
                    if len(seq) - 1 == ranges[1][-1] and \
                                    ranges[0][0] == 0:
                        features[f] = (ranges[1][0], ranges[0][-1])
                    else:
                        raise HydraSequenceException('Feature has more than one location. Could not circularize feature.')
                else:
                    features[f] = ranges
            elif len(ranges) > 2:
                raise HydraSequenceException('Feature has more than one location')
        features = dict(features)
        return features

    def get_feature(self, feature):
        return self.features()[feature]

    def _filter_features(self, fnc):
       features = self.features()
       found_features = set()
       for f in features:
           if fnc(f):
               found_features.add(f)
       return found_features

    def search_features(self, name):
        return self._filter_features(lambda x: x.name==name)

    def _inbounds(self, num):
        if isinstance(num, int):
            num = [num]
        for n in num:
            if n < 0 or n > len(self.seq()) - 1:
                raise HydraSequenceException("Index {} out of bounds ({}, {})".format(n, 0, len(self)))

    def _inbounds_se(self, start, end):
        self._inbounds((start, end))
        if start > end:
            if not self.is_circular():
                raise HydraSequenceException("Start index {} must be less than end index {} for linear sequences.".format(start, end))

    def _ranges(self, lst):
        pos = (j - i for i, j in enumerate(lst))
        t = 0
        for i, els in groupby(pos):
            l = len(list(els))
            el = lst[t]
            t += l
            yield range(el, el + l)

    def splice(self, start, end):
        # TODO: splice the feature as well
        self._inbounds_se(start, end)
        seq_copy = deepcopy(self)
        spliced = seq_copy.get(start, end)
        spliced_start = spliced[0]
        spliced_end = spliced[-1]

        spliced_sequence = Sequence(first_nucleotide=spliced_start)
        cut_features_1 = spliced_start.features
        cut_features_2 = spliced_end.features

        spliced_start.cut_prev()
        spliced_end.cut_next()
        return spliced_sequence

    def get_feature_from_nt(self, nt, feature):
        return nt.rev(stop_criteria=lambda x: feature not in x.features)[::-1][:-1] + nt.fwd( \
            stop_criteria=lambda x: feature not in x.features)

    def cut(self, location):
        pass

    def insert(self, location):
        pass

    def delete(self, location):
        pass

    def get(self, start, end):
        self._inbounds_se(start, end)
        return self.seq(start_nt=self.seq()[start], end_nt=self.seq()[end])

    def __add__(self, other):
        if self.is_circular():
            raise HydraSequenceException("Cannot concatentate circular plasmids. Use .insert instead.")
        new_seq = deepcopy(self)
        other_seq = deepcopy(other)
        new_seq.first_base.find_last().set_next(other_seq[0])
        return new_seq

    def __str__(self):
        seq = self.seq()
        return ''.join([str(x) for x in seq])

    def __repr__(self):
        return '<HydraSequence | Sequence: {}'.format(
            str(self)
        )

    # def __getitem__(self, key):
    #     seq = self.seq()
    #     return seq[key]

    def __len__(self):
        return len(self.seq())
