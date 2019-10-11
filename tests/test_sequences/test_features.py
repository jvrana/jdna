from copy import copy

import pytest

from jdna.sequence import Feature
from jdna.sequence import Nucleotide
from jdna.sequence import Sequence


class TestFeatureProperties:
    @pytest.fixture(scope="function")
    def seq(self, test_str):
        return Sequence(test_str)

    @pytest.fixture(scope="function")
    def basic_feature(self, test_seq):
        f = test_seq.annotate(1, len(test_seq) - 1, "feature", "type")
        return f

    def test_(self, test_seq):
        pass

    def test_sequence_has_features(self, test_seq):
        print("Head: {}".format(test_seq._head))
        len(test_seq)
        test_seq.annotate(1, len(test_seq) - 1, "feature", "type")
        test_seq.annotate(1, len(test_seq) - 1, "feature", "type")
        test_seq.annotate(1, len(test_seq) - 1, "feature", "type")
        assert len(test_seq.features_list) == 3

    def test_sequence_has_feature_positions(self, test_seq, test_str):
        f1 = test_seq.annotate(1, len(test_seq) - 1, "feature", "type")
        f2 = test_seq.annotate(2, len(test_seq) - 2, "feature", "type")
        f3 = test_seq.annotate(3, len(test_seq) - 3, "feature", "type")
        fpos = test_seq.features()
        assert len(fpos) == 3
        assert fpos[f1] == [[1, len(test_str) - 1]]
        assert fpos[f2] == [[2, len(test_str) - 2]]
        assert fpos[f3] == [[3, len(test_str) - 3]]

    # def test_feature_has_nodes(self, basic_feature):
    #     print(basic_feature.nodes)

    # def test_feature_has_segments(self, basic_feature):
    #     print(basic_feature.segments)

    def test_basic_feature_name(self, basic_feature):
        assert basic_feature.name == "feature"

    def test_basic_feature_type(self, basic_feature):
        assert basic_feature.type == "type"

    # def test_basic_feature_start(self, basic_feature):
    #     assert basic_feature.start == 1
    #
    # def test_basic_feature_end(self, basic_feature):
    #     assert basic_feature.end == 25 - 1


class TestFeatureManipulations:
    def test_features_linearizing(self, test_seq):
        test_seq.annotate(0, len(test_seq) - 1, "feature", "type")
        test_seq.circularize()
        test_seq.linearize()

    @pytest.mark.parametrize("c,s,e", [(5, 0, 10)])
    def test_feature_cutting(self, test_seq, c, s, e):
        """
        e.g.
        012345678912345  INDICES
        XXXXXXXXXXX----  FEATURE LOCATION

             V <--CUT HERE
        01234 5678912345  OLD INDEX
        XXXXX XXXXXX----  FEATURE LOCATION
        01234 0123456789  NEW INDEX
        """
        f = Feature(name="feature name", type="feature type")
        test_seq.add_feature(s, e, f)
        assert len(test_seq.features_list) == 1
        fragments = test_seq.cut(c)
        assert len(fragments[0].features_list) == 1
        assert len(fragments[1].features_list) == 1

        feature = list(fragments[0].features_list)[0]
        assert dict(fragments[0].features()) == {feature: [[s, c - 1]]}
        assert dict(fragments[1].features()) == {feature: [[0, e - c]]}

        # print(feature.segments())
        # assert len(feature.segments()) == 2

    def test_feature_fusion(self, test_seq):
        f1 = Feature(name="myfeature", type="myfeature")
        f2 = Feature(name="myfeature", type="myfeature")

        i = 0
        j = 4
        k = 5
        l = 7

        test_seq.add_feature(i, j, f1)
        test_seq.add_feature(k, l, f2)
        assert len(test_seq.features_list) == 2
        Nucleotide.fuse_features(test_seq.get(j), test_seq.get(k))
        assert len(test_seq.features_list) == 1
        assert test_seq.features() == {f1: [[i, l]]}

    def test_feature_fusion_unsuccessful(self, test_seq):
        f1 = Feature(name="myfeature", type="myfeature")
        f2 = Feature(name="myfeature2", type="myfeature")

        i = 0
        j = 4
        k = 5
        l = 7

        test_seq.add_feature(i, j, f1)
        test_seq.add_feature(k, l, f2)
        assert len(test_seq.features_list) == 2
        Nucleotide.fuse_features(test_seq.get(j), test_seq.get(k))
        assert len(test_seq.features_list) == 2

    def test_feature_fusion_cyclic(self, test_seq):
        f1 = Feature(name="name")

        i, j, k, l = len(test_seq) - 5, len(test_seq) - 1, 0, 5

        test_seq.add_feature(i, j, f1)
        test_seq.add_feature(k, l, f1)
        assert test_seq.features() == {f1: [[k, l], [i, j]]}

        test_seq.cyclic = True
        assert test_seq.features() == {f1: [[i, l]]}

        # for start in range(1, len(test_seq)):
        #     for end in range(start, len(test_seq)):
        #         seq_copy = copy(seq)
        #         seq_copy.add_feature(start, end, f1)
        #         for cut in range(start, end):
        #             fragments = seq_copy.cut(cut)
        #             fused = fragments[0].insert(fragments[1], len(fragments[0]))
        #             assert len(fused.feature_positions()) == 1
        #             feature = list(fused.feature_positions().keys())[0]
        #             assert fused.feature_positions()[feature] == [(start, end), (0, end - start)]

    def test_add_feature_to_cyclic(self, test_seq):
        f = Feature(name="feature fail", type="failure")
        with pytest.raises(IndexError):
            test_seq.add_feature(-1, 10, f)
        with pytest.raises(IndexError):
            test_seq.add_feature(10, len(test_seq), f)
        with pytest.raises(IndexError):
            test_seq.add_feature(10, 1, f)
        test_seq.circularize()
        start = 10
        end = 1
        test_seq.add_feature(start, end, f)
        assert test_seq.features() == {f: [[10, 1]]}

    def test_copying_features(self, test_seq):
        f = Feature(name="alfjl", type="lkdjflkj")
        test_seq.add_feature(0, 5, f)
        seq_copy = copy(test_seq)

        assert id(f) not in [id(_f) for _f in list(seq_copy.features().keys())]
        assert len(seq_copy.features()) == 1
        assert f in test_seq.features()
        print(seq_copy.features_list[0])
        print(seq_copy.features())

    def test_multiple_features(self, test_seq):
        f1 = Feature(name="feature 1", type="test")
        f2 = Feature(name="feature 2", type="test")
        s1, e1 = (0, 10)
        s2, e2 = (5, 15)
        test_seq.add_feature(s1, e1, f1)
        test_seq.add_feature(s2, e2, f2)
        assert test_seq.features() == {f1: [[s1, e1]], f2: [[s2, e2]]}
        assert f1 in test_seq.nodes[7].features_list
        assert f2 in test_seq.nodes[7].features_list

    def test_maintain_features_after_reverse(self, test_seq, test_str):
        f1 = Feature(name="feature 1", type="test")
        test_seq.add_feature(5, 10, f1)
        assert len(test_seq.features_list) == 1
        test_seq.reverse()
        # assert str(test_seq) == test_str[::-1]
        assert len(test_seq) == len(test_str)
        print(test_seq.features_list)
        assert (
            len(test_seq.features_list) == 1
        ), "Sequence should have one feature after reversal"
        # assert test_seq.feature_positions() == {f1: [(len(test_seq) - 1 - 10, len(test_seq) - 1 - 5), (5, 0)]}
