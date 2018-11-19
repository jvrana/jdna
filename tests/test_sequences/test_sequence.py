from copy import copy

import pytest

from jdna import Feature, Sequence, Nucleotide


# <editor-fold desc="Basic Tests">
def test_Sequence():
    seq = 'This is a test string for the linked data set.'
    l = Sequence(sequence=seq)
    assert str(l) == seq
    assert l.head.data, 'T'
    assert l.head.find_last().data, '.'
    assert len(l) == len(seq)


def test_cyclic_vs_liner():
    seq = 'This is a test string for the linked data set.'
    l = Sequence(sequence=seq)
    assert ~l.cyclic

    l.circularize()
    assert l.cyclic

    assert str(l) == seq
    assert l.head.data, 'T'
    assert l.head.find_last().data, '.'
    assert len(l) == len(seq)

    l.linearize(1)
    assert ~l.cyclic


def test_reindexing():
    seq = 'This is a test string for the linked data set.'
    l = Sequence(sequence=seq)
    l.circularize()

    new_index = 10

    l.linearize(new_index)
    nodes = l.nodes
    assert str(l) == seq[new_index:] + seq[:new_index]
    assert ~l.cyclic


def test_copy():
    seq = 'This is a test string for the linked data set.'
    l = Sequence(sequence=seq)
    l.name = 'sequence name'
    l_copy = copy(l)
    assert str(l) == str(l_copy)
    assert l.name == l_copy.name
    assert l_copy.nodes[0] is not l.nodes[0]
    for n1, n2 in zip(l.nodes, l_copy.nodes):
        assert n1 is not n2
    assert l.__class__ == l_copy.__class__
    from copy import deepcopy
    with pytest.raises(NotImplementedError):
        deepcopy(l)


def test_insertion():
    seq = 'This is a test string for the linked data set.'
    l = Sequence(sequence=seq)
    for i in range(len(l) + 1):
        l_copy = copy(l)
        insertion = Sequence(sequence='XYZ')
        l_copy.insert(insertion, i)
        assert str(l_copy) == seq[:i] + str(insertion) + seq[i:]


def test_insertion2():
    seq = Sequence(sequence='agtcgagggcagatag')
    seq2 = Sequence(sequence='GTAGGCAG')

    for i in range(len(seq)):
        assert str(copy(seq).insert(copy(seq2), i)) == str(seq)[:i] + str(seq2) + str(seq)[i:]


def test_insertion_index_error():
    seq = 'This is a test string for the linked data set.'
    l = Sequence(sequence=seq)
    with pytest.raises(IndexError):
        l.insert(Sequence(sequence='XYZ'), -1)

    with pytest.raises(IndexError):
        l.insert(Sequence(sequence='XYZ'), len(l) + 1)


def test_circular_cutting():
    seq = 'This is a test string for the linked data set.'
    l = Sequence(sequence=seq)
    l.circularize()
    for cs in [[10, 15, 20], [4, 10, 15]]:
        fragments = [str(x) for x in l.cut(cs)]
        expected = [seq[cs[-1]:] + seq[:cs[0]]]
        expected += [seq[cs[i]:cs[i + 1]] for i in range(len(cs) - 1)]
        assert set(expected) == set(fragments)


def test_linear_cutting():
    seq = 'This is a test string for the linked data set.'
    l = Sequence(sequence=seq)
    for cs in [[10, 15, 20], [4, 10, 15]]:
        fragments = [str(x) for x in l.cut(cs)]
        expected = [seq[:cs[0]]]
        expected += [seq[cs[i]:cs[i + 1]] for i in range(len(cs) - 1)]
        expected += [seq[cs[-1]:]]
        expected = set(expected)
        fragments = set(fragments)
        assert expected == fragments


def test_search():
    query_seq = 'ABCDEFGHIJKLMNOP'
    template = Sequence(sequence=query_seq)
    for i in range(len(template) + 1):
        for j in range(i + 1, len(template) + 1):
            query = Sequence(sequence=query_seq[i:j])
            assert (i, template.nodes[i]) == template.search_all(Sequence(sequence=query_seq[i:j]))[0]

    query = Sequence(sequence='ABCDFG')
    assert [] == template.search_all(query)

    template = Sequence(sequence='XXXXGHHHXHGG')
    query = Sequence(sequence='XX')
    assert 3 == len(template.search_all(query))

    template = Sequence(sequence='longer sentence now for this sentence to find the sentence sentence')
    query = Sequence(sequence='sentence')
    assert 4 == len(template.search_all(query))


def test_circular_search():
    template = Sequence(sequence='XXXXGHHHXHGG')
    template.circularize()
    query = Sequence(sequence='HGGXX')
    assert 1 == len(template.search_all(query))


def test_reverse():
    seq = 'XXXXGHHHXHGG'
    l = Sequence(sequence=seq)
    l.reverse()
    assert str(l) == seq[::-1]

    seq = 'XXXXGHHHXHGG'
    l = Sequence(sequence=seq)
    l.circularize()
    l.reverse()
    assert str(l) == seq[::-1]

    seq = 'XXXXGHHHXHGG'
    l = Sequence(sequence=seq)
    l.reverse()
    l.circularize()
    assert str(l) == seq[::-1]


# </editor-fold>


"""
feature cutting
copying features
replacing features
copy features from

"""


# <editor-fold desc="DNA Manipulation methods">
def test_cutting():
    seq_str = 'agtcgtatgctgcggcgattctgatgctgatgctgatgtcggta'
    seq = Sequence(sequence=seq_str)

    def test_cut(c):
        expected = [seq_str[:c[0]]]
        prev = c[0]
        for loc in c[1:]:
            expected.append(seq_str[prev:loc])
            prev = loc
        expected += [seq_str[c[-1]:]]
        while '' in expected:
            expected.remove('')
        fragments = seq.cut(c)
        assert set([str(x) for x in fragments]) == set(expected)

    test_cut([1, 5, 7])
    test_cut((5, len(seq) - 1))
    test_cut([0, 1, 3, len(seq) - 1, len(seq)])
    test_cut([0, 5, 7])
    test_cut([1, 8, len(seq)])
    assert str(seq) == seq_str


def test_cutting_next():
    seq_str = 'agtcgtatgctgcggcgattctgatgctgatgctgatgtcggta'
    seq = Sequence(sequence=seq_str)

    def test_cut(c):
        expected = [seq_str[:c[0] + 1]]
        prev = c[0]
        for loc in c[1:]:
            expected.append(seq_str[prev + 1:loc + 1])
            prev = loc
        expected += [seq_str[c[-1] + 1:]]
        while '' in expected:
            expected.remove('')
        fragments = seq.cut(c, cut_prev=False)
        assert [str(x) for x in fragments] == expected

    test_cut([1, 5, 7])
    test_cut((5, len(seq) - 1))
    with pytest.raises(IndexError):
        test_cut([0, 1, 3, len(seq) - 1, len(seq)])
    test_cut([0, 5, 7])
    test_cut([1, 8, len(seq) - 1])
    assert str(seq) == seq_str


def test_cyclic_cutting():
    seq_str = 'agtcgtatgctgcggcgattctgatgctgatgctgatgtcggta'
    seq = Sequence(sequence=seq_str)
    seq.circularize()

    def test_cut(c):
        expected = [seq_str[c[-1]:] + seq_str[:c[0]]]
        prev = c[0]
        for loc in c[1:]:
            expected.append(seq_str[prev:loc])
            prev = loc
        while '' in expected:
            expected.remove('')
        fragments = seq.cut(c)
        for f in fragments:
            assert ~f.cyclic
        assert set([str(x) for x in fragments]) == set(expected)

    test_cut([1])
    test_cut([1, 5, 7])
    test_cut((5, len(seq) - 1))
    test_cut([0, 1, 3, len(seq) - 1, len(seq)])
    test_cut([0, 5, 7])
    test_cut([1, 8, len(seq)])
    assert str(seq) == seq_str


def test_complement():
    seq_str = 'AGTCAGTC'
    seq = Sequence(sequence=seq_str)
    seq.complement()
    assert str(seq) == 'TCAGTCAG'
    seq.complement()
    assert str(seq) == seq_str


def test_cyclic_complement():
    seq_str = 'AGTcaGTC'
    seq = Sequence(sequence=seq_str)
    seq.circularize()
    seq.complement()
    assert str(seq) == 'TCAgtCAG'


def test_reverse_complement():
    seq_str = 'AGTCAGTC'
    seq = Sequence(sequence=seq_str)
    seq.reverse_complement()
    assert str(seq) == 'TCAGTCAG'[::-1]
    seq.reverse_complement()
    assert str(seq) == seq_str


def test_cyclic_reverse_complement():
    seq_str = 'AGTCAGTC'
    seq = Sequence(sequence=seq_str)
    seq.circularize()
    seq.reverse_complement()
    assert str(seq) == 'TCAGTCAG'[::-1]


# </editor-fold>

def test_chop():
    seq = Sequence(sequence='agtcgtatgctgcggcgattctgatgctgatgctgatgtcggta')
    for i in range(len(seq)):
        assert str(seq.chop_off_fiveprime(i)) == str(seq)[i:]

    for i in range(len(seq)):
        assert str(seq.chop_off_threeprime(i)) == str(seq)[:i + 1]


class TestFeature(object):
    class TestBasicFeature(object):

        @pytest.fixture(scope='function')
        def seq(self):
            return Sequence(sequence='abcdefghijklmnop')

        @pytest.fixture(scope="function")
        def basic_feature(self, seq):
            f = seq.create_feature('feature', 'type', 1, len(seq) - 1)
            return f

        def test_(self, seq):
            pass

        def test_sequence_has_features(self, seq):
            print("Head: {}".format(seq._head))
            len(seq)
            seq.create_feature('feature', 'type', 1, len(seq) - 1)
            seq.create_feature('feature', 'type', 1, len(seq) - 1)
            seq.create_feature('feature', 'type', 1, len(seq) - 1)
            assert len(seq.features) == 3

        def test_sequence_has_feature_positions(self, seq):
            f1 = seq.create_feature('feature', 'type', 1, len(seq) - 1)
            f2 = seq.create_feature('feature', 'type', 2, len(seq) - 2)
            f3 = seq.create_feature('feature', 'type', 3, len(seq) - 3)
            fpos = seq.feature_positions()
            assert len(fpos) == 3
            assert fpos[f1] == [[1, 15]]
            assert fpos[f2] == [[2, 14]]
            assert fpos[f3] == [[3, 13]]

        def test_feature_has_nodes(self, basic_feature):
            print(basic_feature.nodes)

        def test_feature_has_segments(self, basic_feature):
            print(basic_feature.segments)

        def test_basic_feature_name(self, basic_feature):
            assert basic_feature.name == 'feature'

        def test_basic_feature_type(self, basic_feature):
            assert basic_feature.type == 'type'

        def test_basic_feature_start(self, basic_feature):
            assert basic_feature.start == 1

        def test_basic_feature_end(self, basic_feature):
            assert basic_feature.end == 25 - 1

    def test_features_linearizing(self):
        seq = Sequence(sequence='agttggagcg')
        seq.create_feature('feature', 'type', 0, len(seq) - 1)
        seq.circularize()
        seq.linearize()

    def test_feature_cutting(self):
        """

        012345678912345  INDICES
        XXXXXXXXXXX----  FEATURE LOCATION

             V <--CUT HERE
        01234 5678912345  OLD INDEX
        XXXXX XXXXXX----  FEATURE LOCATION
        01234 0123456789  NEW INDEX
        """
        seq = Sequence(sequence='agtcgtatgctgcggcgattctgatgctgatgctgatgtcggta')
        f = Feature(name='feature name', type='feature type')
        s, e = 0, 10
        seq.add_feature(s, e, f)
        assert len(seq.features) == 1
        c = 5
        fragments = seq.cut(c)
        assert len(fragments[0].features) == 1
        assert len(fragments[1].features) == 1

        feature = list(fragments[0].features)[0]
        assert dict(fragments[0].feature_positions()) == {feature: [[s, c - 1]]}
        assert dict(fragments[1].feature_positions()) == {feature: [[0, e - c]]}

    def test_feature_fusion(self):
        seq = Sequence(sequence='abcdefghijkl')
        f1 = Feature(name='myfeature', type='myfeature')
        f2 = Feature(name='myfeature', type='myfeature')

        i = 0
        j = 4
        k = 5
        l = 7

        seq.add_feature(i, j, f1)
        seq.add_feature(k, l, f2)
        assert len(seq.features) == 2
        Nucleotide.fuse_features(seq.get(j), seq.get(k))
        assert len(seq.features) == 1
        assert seq.feature_positions() == {
            f1: [[i, l]]
        }

    def test_feature_fusion_unsuccessful(self):
        seq = Sequence(sequence='abcdefghijkl')
        f1 = Feature(name='myfeature', type='myfeature')
        f2 = Feature(name='myfeature2', type='myfeature')

        i = 0
        j = 4
        k = 5
        l = 7

        seq.add_feature(i, j, f1)
        seq.add_feature(k, l, f2)
        assert len(seq.features) == 2
        Nucleotide.fuse_features(seq.get(j), seq.get(k))
        assert len(seq.features) == 2

    def test_feature_fusion_by_setting(self):
        seq1 = Sequence(sequence='a'*10)
        seq2 = Sequence(sequence='a'*10)

        f1 = Feature(name='myfeature')
        f2 = Feature(name='myfeature')

        seq1.add_feature(5, 9, f1)
        seq2.add_feature(0, 5, f2)

        product = seq1.fuse(seq2)
        assert len(product.features) == 1
        assert product.feature_positions() == {
            f1: [[5, 15]]
        }

    def test_feature_fusion_by_setting(self):
        seq1 = Sequence(sequence='a'*10)
        seq2 = Sequence(sequence='a'*10)

        f1 = Feature(name='myfeature')
        f2 = Feature(name='myfeature')

        seq1.add_feature(5, 9, f1)
        seq2.add_feature(1, 5, f2)

        product = seq1.fuse(seq2)
        assert len(product.features) == 2

    def test_feature_fusion_cyclic(self):
        seq1 = Sequence(sequence='a'*20)
        f1 = Feature(name='name')
        seq1.add_feature(10, 19, f1)
        seq1.add_feature(0, 5, f1)
        assert seq1.feature_positions() == {
            f1: [[0, 5], [10, 19]]
        }

        seq1.cyclic = True
        assert seq1.feature_positions() == {
            f1: [[10, 5]]
        }


        # for start in range(1, len(seq)):
        #     for end in range(start, len(seq)):
        #         seq_copy = copy(seq)
        #         seq_copy.add_feature(start, end, f1)
        #         for cut in range(start, end):
        #             fragments = seq_copy.cut(cut)
        #             fused = fragments[0].insert(fragments[1], len(fragments[0]))
        #             assert len(fused.feature_positions()) == 1
        #             feature = list(fused.feature_positions().keys())[0]
        #             assert fused.feature_positions()[feature] == [(start, end), (0, end - start)]

    def test_add_feature_to_cyclic(self):
        seq = Sequence(sequence='agtcgtatgctgcggcgattctgatgctgatgctgatgtcggta')
        f = Feature(name='feature fail', type='failure')
        with pytest.raises(IndexError):
            seq.add_feature(-1, 10, f)
        with pytest.raises(IndexError):
            seq.add_feature(10, len(seq), f)
        with pytest.raises(IndexError):
            seq.add_feature(10, 1, f)
        seq.circularize()
        start = 10
        end = 1
        seq.add_feature(start, end, f)
        assert seq.feature_positions() == {f: [[10, 1]]}

    def test_add_and_feature_positions(self):
        def s(seq_str):
            return Sequence(sequence=seq_str)

        seq = s('agtcgtatgctgcggcgattctgatgctgatgctgatgtcggta')
        for j in range(0, len(seq)):
            for k in range(j, len(seq), 3):
                f = Feature(name='new feature', type='feature type')
                seq_copy = copy(seq)
                for l in seq_copy.nodes:
                    assert l.features == set()
                assert seq_copy.feature_positions() == {}
                seq_copy.add_feature(j, k, f)
                for i in range(j, k + 1):
                    assert f in seq_copy.nodes[i].features

                features = seq_copy.feature_positions()
                assert len(features) == 1
                assert features == {f: [(j, k), (0, k - j)]}

    def test_copying_features(self):
        seq = Sequence(sequence="AGTCAGTCGGA")
        f = Feature(name='alfjl', type='lkdjflkj')
        seq.add_feature(0, 5, f)
        seq_copy = copy(seq)
        assert f not in list(seq_copy.feature_positions().keys())
        assert len(seq_copy.feature_positions()) == 1
        assert f in seq.feature_positions()

    def test_multiple_features(self):
        def s(seq_str):
            return Sequence(sequence=seq_str)

        seq = s('agtcgtatgctgcggcgattctgatgctgatgctgatgtcggta')
        f1 = Feature(name='feature 1', type='test')
        f2 = Feature(name='feature 2', type='test')
        s1, e1 = (0, 10)
        s2, e2 = (5, 15)
        seq.add_feature(s1, e1, f1)
        seq.add_feature(s2, e2, f2)
        assert seq.feature_positions() == {f1: [(s1, e1), (0, e1 - s1)], f2: [(s2, e2), (0, e2 - s2)]}
        assert f1 in seq.nodes[7].features
        assert f2 in seq.nodes[7].features

    def test_reindexed_features(self):
        def s(seq_str):
            return Sequence(sequence=seq_str)

        seq = s('agtcgtatgctgcggcgattctgatgctgatgctgatgtcggta')
        f1 = Feature(name='feature 1', type='test')
        f2 = Feature(name='feature 2', type='test')
        s1, e1 = (0, 10)
        s2, e2 = (5, 15)
        si1, si2 = 2, 3
        seq.add_feature(s1, e1, f1, start_index=si1)
        seq.add_feature(s2, e2, f2, start_index=si2)
        assert seq.feature_positions() == {f1: [(s1, e1), (si1, si1 + e1 - s1)], f2: [(s2, e2), (si2, si2 + e2 - s2)]}

    def test_maintain_features_after_reverse(self):
        seq = Sequence(sequence='agtcgtatgctgcggcgattctgatgctgatgctgatgtcggta')
        f1 = Feature(name='feature 1', type='test')
        seq.add_feature(5, 10, f1)
        seq.reverse()
        assert seq.feature_positions() == {f1: [(len(seq) - 1 - 10, len(seq) - 1 - 5), (5, 0)]}
