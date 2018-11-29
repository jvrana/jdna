import pytest

from jdna.reaction import Reaction
from jdna.sequence import Sequence


@pytest.fixture(scope='function')
def template():
    return Sequence.random(500)


@pytest.mark.parametrize('o1', [10, 0])
@pytest.mark.parametrize('o2', [10, 0])
@pytest.mark.parametrize('primer1_inset', [
    pytest.param(30, id="inset1=30bp"),
    pytest.param(15, id="inset1=15bp"),
])
@pytest.mark.parametrize('primer2_inset', [
    pytest.param(30, id="inset2=30bp"),
    pytest.param(15, id="inset2=15bp"),
])
@pytest.mark.parametrize('circular_template', [False, True])
def test_pcr(primer1_inset, primer2_inset, o1, o2, circular_template):
    template = Sequence.random(100)

    primer_len = 20
    p1_start = primer1_inset
    p1_end = primer1_inset + primer_len

    p2_end = len(template) - primer2_inset
    p2_start = p2_end - primer_len

    print(p1_start)
    print(p1_end)
    print(p2_start)
    print(p2_end)

    p1 = Sequence('N'*o1) + template[p1_start:p1_end]
    p2 = Sequence('N'*o2) + template[p2_start:p2_end].reverse_complement()

    template.circularize()
    if circular_template:
        template.reindex(int((len(template) - primer2_inset - primer1_inset)/2))

    products = Reaction.pcr(template, p1, p2)
    assert len(products) == 1
    expected_len = len(template) - primer1_inset - primer2_inset + o1 + o2
    assert len(products[0]) == expected_len, "Product should be {} long".format(expected_len)

    # p1 = Sequence('N'*o1) + template[pos1:pos1+l1]
    # p2 = Sequence('N'*o2) + template[pos2-l2:pos2+1].reverse_complement()
    #
    # products = Reaction.pcr(template, p1, p2)
    # assert len(products) == 1
    # expected_len = o1 + o2 + (pos2-pos1)
    # assert len(products[0]) == expected_len
    # assert str(products[0]) == 'N'*o1 + str(template)[pos1:pos2] + 'N'*o2


@pytest.fixture(scope='function')
def seq():
    return Sequence.random(300)

@pytest.fixture(scope='function')
def generate_sequences(seq):
    """
    A number of overlapping sequences
    :param test_str:
    :type test_str:
    :return:
    :rtype:
    """
    def generate_sequences(num_fragments=3, overhang=20, cyclic=True):
        sequences = []
        indices = list(range(0, len(seq), int(len(seq)/(num_fragments))))
        j = 0
        for i, j in zip(indices[:-1], indices[1:]):
            sequences.append(seq[i:j+overhang])
        sequences.append(seq[j:])
        if cyclic:
            sequences[-1].fuse_in_place(seq[:overhang])
        for i, s in enumerate(sequences):
            s.name = str(i)
        return sequences
    return generate_sequences


def test_(generate_sequences):
    seqs = generate_sequences(3, cyclic=False)
    bindings = Reaction.anneal_sequences(seqs)
    for b in bindings:
        print()
        print(b.primer.global_id)
        print(b.template.global_id)
        print(b.position)


def test_anneal():
    template = Sequence.random(200)

    anneal = template[50:80]
    overhang = Sequence('N'*20)
    primer = overhang + anneal

    bindings = template.anneal(primer)
    for b in bindings:
        print(b)


@pytest.mark.parametrize('cyclic', ['cyclic', False])
@pytest.mark.parametrize('overhang_length', [20, 30, 15, 10])
@pytest.mark.parametrize('num_fragments', [1, 3])
@pytest.mark.parametrize('reverse_first', [False, 'reverse_first'])
def test_anneal_sequences(seq, num_fragments, generate_sequences, cyclic, overhang_length, reverse_first):
    cyclic = cyclic == 'cyclic'
    reverse_first = reverse_first == 'reverse_first'
    sequences = generate_sequences(num_fragments=num_fragments, cyclic=cyclic, overhang=overhang_length)

    if reverse_first:
        sequences[0].reverse_complement()
    bindings = list(Reaction.anneal_sequences(sequences))

    if overhang_length >= Sequence.DEFAULTS.MIN_ANNEAL_BASES:
        assert len(bindings) == (num_fragments - int(not cyclic)) * 2
    else:
        assert len(bindings) == 0

    for b in bindings:
        assert b.position.span[1] - b.position.span[0] + 1 == overhang_length, "Length of overhang should be {}".format(overhang_length)

    for i, b in enumerate(bindings):
        if reverse_first:
            if (cyclic and i in [0, 1]) or (not cyclic and i == 0):
                assert str(b.position.template_anneal.reverse_complement()) in str(seq), \
                'Reverse sequence ({}) {} not in sequence\n{}'.format(i, b.position.template_anneal.reverse_complement(), seq)
        else:
            assert str(b.position.template_anneal) in str(seq), 'Sequence ({}) {} not in template\n{}'.format(i, b.position.template_anneal, seq)


@pytest.mark.parametrize('cyclic', [pytest.param(False, id='linear'), pytest.param(True, id='circular'),])
@pytest.mark.parametrize('try_reverse_complement', [False, True])
@pytest.mark.parametrize('bind_reverse_complement', [False, True])
def test_interaction_graph(cyclic, seq, generate_sequences, try_reverse_complement, bind_reverse_complement):
    num_fragments = 4
    sequences = generate_sequences(num_fragments, cyclic=cyclic)
    if try_reverse_complement == 'try_reverse_complement':
        sequences[1].reverse_complement()
    G = Reaction.interaction_graph(sequences, bind_reverse_complement=bind_reverse_complement)

    edges = list(G.edges)
    print(len(edges))
    if bind_reverse_complement:
        if cyclic:
            assert len(edges) == 4 * 2
        else:
            assert len(edges) == 3 * 2
    else:
        if cyclic:
            assert len(edges) == 4
        else:
            assert len(edges) == 3

@pytest.mark.parametrize('cyclic', [pytest.param(False, id='linear'), pytest.param(True, id='circular'),])
def test_linear_paths(generate_sequences, cyclic):
    sequences = generate_sequences(4, cyclic=cyclic)
    G = Reaction.interaction_graph(sequences, bind_reverse_complement=True)
    paths = Reaction.linear_paths(G)

    if not cyclic:
        assert len(paths) == 2
        for p in paths:
            assert len(p) == 4
    else:
        assert len(paths) == 0

@pytest.mark.parametrize('cyclic', [pytest.param(False, id='linear'), pytest.param(True, id='circular'),])
def test_cyclic_paths(generate_sequences, cyclic):
    sequences = generate_sequences(4, cyclic=cyclic)
    sequences[0].reverse_complement()
    G = Reaction.interaction_graph(sequences, bind_reverse_complement=True)
    paths = Reaction.cyclic_paths(G)
    for p in paths:
        print(p)
    for e in G.edges:
        print(e)
    if cyclic:
        assert len(paths) == 2
        for p in paths:
            assert len(p) == 4
    else:
        assert len(paths) == 0

@pytest.mark.parametrize('reverse_complement', [
    [],
    [0],
    [1],
    [0, 1],
    [0, 2],
])
def test_linear_assemblies(seq, generate_sequences, reverse_complement):
    sequences = generate_sequences(4, cyclic=False)
    for rc in reverse_complement:
        sequences[rc].reverse_complement()
    assemblies = Reaction.linear_assemblies(sequences)

    assert len(assemblies) == 2

    for a in assemblies:
        p = a.product
        assert not p.cyclic
        assert len(p) == len(seq)
        assert str(p) == str(seq) or str(p.copy().reverse_complement()) == str(seq)

@pytest.mark.parametrize('reverse_complement', [
    [],
    [0],
    [1],
    [2],
    [0, 1],
    [0, 2],
])
def test_cyclic_assemblies(seq, generate_sequences, reverse_complement):
    sequences = generate_sequences(4, cyclic=True)
    for rc in reverse_complement:
        sequences[rc].reverse_complement()

    assemblies = Reaction.cyclic_assemblies(sequences, depth=50)

    assert len(assemblies) == 2
    for a in assemblies:
        p = a.product
        assert p.cyclic
        expected = seq.copy().circularize()
        rc_expected = expected.copy().reverse_complement()
        g = [p.compare(expected), p.compare(rc_expected)]
        a.align()
        assert any(g)