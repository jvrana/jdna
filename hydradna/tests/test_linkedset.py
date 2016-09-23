from nose.tools import *

from hydradna import *



def test_linkedset():
    seq = 'This is a test string for the linked data set.'
    l = LinkedSet(data_sequence=seq)
    assert_equal(str(l), seq)
    assert_equal(l.get_first().data, 'T')
    assert_equal(l.get_first().find_last().data, '.')


    # len
    assert_equal(len(l), len(seq))


    # Cyclic vs linear
    assert_false(l.is_cyclic())
    l.make_cyclic()
    assert_true(l.is_cyclic())
    l.linearize(0)
    assert_false(l.is_cyclic())
    assert_equal(str(l), seq)


    # Re-index
    l.make_cyclic()
    l.reindex(0)
    assert_equal(str(l), seq)
    for i in range(len(l)):
        l_copy = deepcopy(l)
        l_copy.reindex(i)
        assert_equal(str(l_copy), seq[i:] + seq[:i])


    # Insertion
    for i in range(len(l)):
        l_copy = deepcopy(l)
        insertion = LinkedSet(data_sequence='XYZ')
        l_copy.insert(insertion, i)
        assert_equal(str(l_copy), seq[:i] + str(insertion) + seq[i:])

    l_copy = deepcopy(l_copy)
    assert_raises(IndexError, l_copy.insert, LinkedSet(data_sequence='XYZ'), len(l_copy)+1)
    assert_raises(IndexError, l_copy.insert, LinkedSet(data_sequence='XYZ'), -1)


    # Cutting
    l_copy = deepcopy(l)
    for cs in [[10, 15, 20], [4, 10, 15]]:
        fragments = [str(x) for x in l_copy.cut(cs)]
        expected = [seq[cs[-1]:] + seq[:cs[0]]]
        expected += [ seq[cs[i]:cs[i+1]] for i in range(len(cs)-1) ]
        for e, f in zip(expected, fragments):
            assert_true(str(e) == str(f))

    assert_raises(IndexError, l_copy.cut, len(l_copy)+1)
    assert_raises(IndexError, l_copy.cut, -1)


    l_copy = deepcopy(l)
    l_copy.linearize(0)
    for cs in [[10, 15, 20], [4, 10, 15]]:
        fragments = [str(x) for x in l_copy.cut(cs)]
        expected = [seq[:cs[0]]]
        expected += [ seq[cs[i]:cs[i+1]] for i in range(len(cs)-1) ]
        expected += [seq[cs[-1]:]]
        for e, f in zip(expected, fragments):
            assert_true(str(e) == str(f))
            assert_equal(len(fragments), len(cs)+1)

    #Search
    i, j = 10, 15
    query_seq = 'ABCDEFGHIJKLMNOP'
    template = LinkedSet(data_sequence=query_seq)
    for i in range(len(template)+1):
        for j in range(i+1, len(template)+1):
            query = LinkedSet(data_sequence=query_seq[i:j])
            assert_equal((i, template.get()[i]), template.search_all(LinkedSet(data_sequence=query_seq[i:j]))[0])

    query = LinkedSet(data_sequence='ABCDFG')
    assert_equal([], template.search_all(query))

    template = LinkedSet(data_sequence='XXXXGHHHXHGG')
    query = LinkedSet(data_sequence='XX')
    assert_equal(3, len(template.search_all(query)))

    template = LinkedSet(data_sequence='longer sentence now for this sentence to find the sentence sentence')
    query = LinkedSet(data_sequence='sentence')
    assert_equal(4, len(template.search_all(query)))


    #Search Circular
    template = LinkedSet(data_sequence='XXXXGHHHXHGG')
    template.make_cyclic()
    query = LinkedSet(data_sequence='HGGXX')
    assert_equal(1, len(template.search_all(query)))



    # Reverse
    l_copy = deepcopy(l)
    l_copy.reverse()
    assert_true(seq[::-1], str(l_copy.reverse()))

    l_copy.make_cyclic()
    l_copy.reverse()

    #template.linearize()