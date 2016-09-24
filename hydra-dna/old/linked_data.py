from collections import defaultdict
from itertools import groupby
from copy import deepcopy
import warnings
import uuid
import itertools

MIN_BASES = 13

class HydraSequenceException(Exception):
    """ Generic exception for HydraSequence errors """


class Link(object):
    """
    A link represents a single piece of data in a linked sequence.
    Each link can be connected to another link through the _next
    and _prev pointers (representing backbone bonds).

    Available methods:
       get_next
       get_prev
       add_next
       add_prev
       cut_next
       cut_prev
       remove_next
       remove_prev
       find_first
       find_last
       fwd
       rev
    """

    def __init__(self, data):
        self.data = data
        self.__next = None
        self.__prev = None

    def next(self):
        return self.__next

    def prev(self):
        return self.__prev

    def add_next(self, data):
        new_link = Link(data)
        self.set_next(new_link)
        return new_link

    def add_prev(self, data):
        new_link = Link(data)
        self.set_prev(new_link)
        return new_link

    def cut_next(self):
        next_link = self.next()
        if next_link is not None:
            next_link.__assign_prev(None)
        self.__assign_next(None)
        return next_link

    def cut_prev(self):
        prev_link = self.prev()
        if prev_link is not None:
            prev_link.__assign_next(None)
        self.__assign_prev(None)
        return prev_link

    def remove(self):
        next_link = self.next()
        prev_link = self.prev()
        if next_link is not None:
            next_link.__assign_prev(prev_link)
        if prev_link is not None:
            prev_link.__assign_next(next_link)
        return

    def swap(self):
        temp = self.__next
        self.__next = self.__prev
        self.__prev = temp

    def __assign_next(self, link):
        self.__next = link

    def __assign_prev(self, link):
        self.__prev = link

    def set_next(self, link):
        if link is not None:
            link.__assign_prev(self)
        self.__assign_next(link)

    def set_prev(self, link):
        if link is not None:
            link.__assign_next(self)
        self.__assign_prev(link)

    def make_cyclic(self):
        if not self.is_cyclic():
            first = self.find_first()
            last = self.find_last()
            last.set_next(first)

    def is_cyclic(self):
        visited = set()
        curr = self
        while curr:
            if curr in visited:
                return True
            visited.add(curr)
            curr = curr.next()
        return False

    def _propogate(self, next_method, stop=None, stop_criteria=None):
        visited = []
        nxt = self
        while nxt is not None:
            if stop_criteria is not None and stop_criteria(nxt):
                break
            visited.append(nxt)
            if nxt is stop:
                break
            nxt = next_method(nxt)
            if nxt is visited[0]:
                break
        return visited

    def fwd(self, stop_link=None, stop_criteria=None):
        return self._propogate(
            lambda x: x.next(),
            stop=stop_link,
            stop_criteria=stop_criteria)

    def rev(self, stop_link=None, stop_criteria=None):
        return self._propogate(
            lambda x: x.prev(),
            stop=stop_link,
            stop_criteria=stop_criteria)

    def find_first(self):
        links = self.rev()
        return links[-1]

    def find_last(self):
        links = self.fwd()
        return links[-1]

    def _longest_match(self, y, next_method):
        x1 = self
        x2 = y
        longest_match = []
        while x1 and x2:
            if x1.equivalent(x2):
                longest_match.append((x1, x2))
                x1 = next_method(x1)
                x2 = next_method(x2)
            else:
                break
        return longest_match

    def _complete_match(self, y, next_method):
        l = self._longest_match(y, next_method)
        if l == []:
            return False
        t1, t2 = self._longest_match(y, next_method)[-1]
        return not (next_method(t1) and next_method(t2))

    def complete_match_fwd(self, y):
        return self._complete_match(y, lambda x: x.next())

    def complete_match_rev(self, y):
        return self._complete_match(y, lambda x: x.prev())

    def equivalent(self, other):
        return self.data == other.data

    def __repr__(self):
        return str(self)

    def __str__(self):
        return str(self.data)


class LinkedSet(object):

    def __init__(self, first=None, data_sequence=None):
        if data_sequence is not None:
            self.initialize(data_sequence)
        elif first is not None:
            self.first = first

    def initialize(self, data_sequence):
        self.first = Link(data_sequence[0])
        current = self.first
        for d in data_sequence[1:]:
            new = Link(d)
            current.set_next(new)
            current = new

    def get_first(self):
        if self.is_cyclic():
            return self.first
        first = self.first.find_first()
        self.first = first
        return self.first

    def make_cyclic(self):
        return self.get_first().make_cyclic()

    def is_cyclic(self):
        visited = set()
        curr = self.first
        while curr:
            if curr in visited:
                return True
            visited.add(curr)
            curr = curr.next()
        return False

    def linearize(self, i=0):
        this_i = self.get()[i]
        this_i.cut_prev()
        return this_i

    def get(self):
        return self.get_first().fwd()

    def cut(self, i, cut_prev=True):
        if isinstance(i, tuple):
            i = list(i)
        if isinstance(i, int):
            i = list(set([i]))
        i = list(set(i))
        i.sort()
        if len(self) in i:
            if cut_prev:
                i.remove(len(self))
        self._inbounds(i)
        self_copy = deepcopy(self)
        all_links = self_copy.get()
        for cut_loc in i:
            link = all_links[cut_loc]
            if cut_prev:
                link.cut_prev()
            else:
                link.cut_next()
        return LinkedSet._group_links(all_links)

    #TODO: Speed up with set
    @staticmethod
    def _group_links(links):
        unique_first_links = []
        for link in links:
            first_link = link.find_first()
            if first_link not in unique_first_links:
                unique_first_links.append(first_link)
        return [LinkedSet(first=l) for l in unique_first_links]

    def insert(self, linkedset, i, copy_insertion=True):
        if i == len(self.get()):
            pass
        else:
            self._inbounds(i)
        if linkedset.is_cyclic():
            raise TypeError("Cannot insert a cyclic sequence")
        if copy_insertion:
            linkedset = deepcopy(linkedset)
        #TODO: This copies the insertion sequence, you want that?
        if i == len(self.get()):
            loc2 = None
            loc1 = self.get()[i-1]
        else:
            loc2 = self.get()[i]
            loc1 = loc2.prev()
        first = linkedset.get()[0]
        last = linkedset.get()[-1]
        first.set_prev(loc1)
        last.set_next(loc2)
        if i == 0: # Special case in which user inserts sequence in front of their sequence; they probably intend to re-index it
            self.first = first
        return self

    def remove(self, i):
        self._inbounds(i)
        return self.get()[i].remove()

    def reindex(self, i):
        self._inbounds(i)
        if not self.is_cyclic():
            raise TypeError("Cannot re-index a linear linked set")
        self.first = self.get()[i]

    def _inbounds(self, num):
        if isinstance(num, int):
            num = [num]
        for n in num:
            mn = 0
            mx = len(self.get()) - 1
            if n < 0 or n > mx:
                raise IndexError("Index {} out of acceptable bounds ({}, {})".format(n, mn, mx))

    def search_all(self, query):
        curr_link = self.get_first()
        q_link = query.get_first()
        i = 0
        found = []
        visited = set()
        while curr_link and curr_link not in visited:
            visited.add(curr_link)
            if curr_link.complete_match_fwd(q_link):
                found.append((i, curr_link))
            curr_link = curr_link.next()
            i += 1
        return found

    def reverse(self, in_place=True):
        x = self
        if not in_place:
            x = deepcopy(x)
        for s in x.get():
            s.swap()
        if x.is_cyclic():
            x.reindex(1)
        return x

    def slice(self, i, j, fwd=True):
        links = self.get()
        start = links[i]
        stop = links[j]
        method = start.fwd
        if not fwd:
            method = start.rev
        sec_links = method(stop_link=stop)
        if sec_links[-1] is stop:
            return sec_links
        else:
            raise IndexError("Improper indices for linkedset.")

    def __reversed__(self):
        l_copy = deepcopy(self)
        for s in l_copy.get():
            s.swap()
        return l_copy

    def __len__(self):
        return len(self.get())

    def __iter__(self):
        current = self.first
        while current is not None:
            yield current
            current = current.next()

    def __repr__(self):
        return str(self)

    def __str__(self):
        return ''.join(str(x) for x in self.get())


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
        next = None
        if cut_prev:
            next = super(Nucleotide, self).cut_prev()
        else:
            next = super(Nucleotide, self).cut_next()

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
                print f1, f2
                if f1.name == f2.name:
                    if f1.end + 1 == f2.start:

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
        return self.data.upper == other.data.upper

class Feature(object):

    def __init__(self, name, type):
        self.name = name
        self.type = type
        self.start = 0
        self.end = None

    def __str__(self):
        return "seq_feature ({}[{}:{}])".format(self.name, self.start, self.end)

    def __repr__(self):
        return str(self)


class Sequence(LinkedSet):

    def __init__(self, first_nt=None, sequence=None):
        super(Sequence, self).__init__(first=first_nt, data_sequence=sequence)

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
        curr = start
        while curr is not end:
            curr.features.add(feature)
            curr = curr.next()
        curr.features.add(feature)

    def get_features(self):
        features = defaultdict(list)
        for i, x in enumerate(self.get()):
            for f in x.features:
                features[f].append(i)
        return features

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

    def complement(self, in_place=True):
        x = self
        if not in_place:
            x = deepcopy(x)
        curr = x.get_first()
        visited = set()
        while curr and curr not in visited:
            visited.add(curr)
            curr.to_complement()
            curr = curr.next()
        return x

    def reverse_complement(self, in_place=True):
        x = self
        if not in_place:
            x = deepcopy(x)
        x.reverse()
        x.complement()
        return x

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

class Reaction(object):

    @staticmethod
    def _anneal(template, primer, min_bases=MIN_BASES, threeprime=True):
        t = template.get_first()
        p = primer.get_first()
        if threeprime:
            t = t.find_last()
            p = p.find_last()
        visited = set()
        i = 0
        matches = []
        while t and t not in visited:
            next_method = lambda x: x.next()
            if threeprime:
                next_method = lambda x: x.prev()
            l = t._longest_match(p, next_method)
            if len(l) >= min_bases:
                matches.append((i, len(l)))
            visited.add(t)
            t = next_method(t)
            i += 1
        return matches

    @staticmethod
    def anneal_fiveprime(template, primer, min_bases=MIN_BASES):
        return Reaction._anneal(template, primer, min_bases=min_bases, threeprime=False)

    @staticmethod
    def anneal_threeprime(template, primer, min_bases=MIN_BASES):
        matches = Reaction._anneal(template, primer, min_bases=min_bases, threeprime=True)
        return [(len(template)-m[0], m[1]) for m in matches]

    @staticmethod
    def anneal_primer(template, primer, min_bases=MIN_BASES):
        rc_template = template.reverse_complement(in_place=False)
        fwd_matches = Reaction.anneal_threeprime(template, primer, min_bases=min_bases)
        # for f in fwd_matches:
        #     primer.cut
        rev_matches = Reaction.anneal_threeprime(rc_template, primer, min_bases=min_bases)

        def ca(matches):
            return [dict(pos=m[0],len=m[1],tm=Reaction.tm(primer.cut(len(primer) - m[1])[-1])) for m in matches]

        anneal = dict(F=ca(fwd_matches), R=ca(rev_matches))

        return anneal

    @staticmethod
    def pcr(template, p1, p2, min_bases=MIN_BASES):
        template = deepcopy(template)

        # Annealing locations for p1
        ann1 = Reaction.anneal_primer(template, p1, min_bases=min_bases)

        # Annealing locations for p2
        ann2 = Reaction.anneal_primer(template, p2, min_bases=min_bases)

        # Group fwd and rev matches
        fwd_matches = ann1['F'] + ann2['F']
        rev_matches = ann1['R'] + ann2['R']

        # All possible paired matches
        pairs = itertools.product(fwd_matches, rev_matches)
        products = []
        for fwd, rev in pairs:
            # Cut the template
            s = fwd['pos']-fwd['len']
            e = len(template)-rev['pos']+rev['len']-1
            start = template.get()[s]
            end = template.get()[e]
            start.cut_prev()
            end.cut_next()
            product = Sequence(first_nt=start)

            # Fuse the overhangs
            overhang1 = p1.cut(len(p1) - fwd['len'])[0]
            overhang2 = p2.cut(len(p2) - rev['len'])[0]
            product.insert(overhang1, 0)
            product.insert(overhang2.reverse_complement(), len(product))
            products.append(product)
        return products

    @staticmethod
    def assembly_cycles(fragments, max_homology, min_homology):
        def complete_match(match):
            return match[0][0] == match[0][1]

        def less_than_max_homology(match):
            return match[0][0] < max_homology

        def one_match(match):
            return len(match) == 1

        def pass_conditions(match):
            return one_match(match) \
                   and less_than_max_homology(match) \
                   and complete_match(match)

        fragments += [deepcopy(f).reverse_complement(in_place=True) for f in fragments[1:]]
        pairs = itertools.permutations(fragments, 2)
        graph = []
        for pair in pairs:
            match = Reaction.anneal_threeprime(*pair, min_bases=min_homology)
            if match and pass_conditions(match):
                graph.append( pair )
        cycles = Utilities.Graph.find_cycles(graph)
        return cycles

    #TODO: make a linear assembly
    @staticmethod
    def cyclic_assembly(fragments, max_homology=60, min_homology=13):
        cycles = Reaction.assembly_cycles(fragments, max_homology=max_homology, min_homology=min_homology)
        if len(cycles) > 1:
            warnings.warn("More than one assembly found.")
        products = []
        for cy in cycles:
            # Copy fragments in cycle
            cy = deepcopy(cy)

            # Pair fragments for assembly
            x1 = cy
            x2 = x1[1:] + x1[:1]
            pairs = zip(x1, x2)

            # Cut 5' ends of fragments according to homology
            for right, left in pairs:
                match = Reaction.anneal_threeprime(right, left)[0]
                nt = right.get()[match[0]]
                nt.cut_prev()
                right.first = nt

            # Fuse all fragments
            for p in pairs:
                right, left = p
                left.fuse(right)

            products.append(x1[0])
        return products

    @staticmethod
    def tm(sequence):
        if not (isinstance(sequence, LinkedSet) or isinstance(sequence, str)):
            raise TypeError("Cannot calculate tm of {} type for {}".format(str(type(sequence)), sequence))
        s = str(sequence).lower()
        gc = s.count('g') + s.count('c')
        return 64.9 + 41 * (gc - 16.4) / len(s)


class Utilities:

    class Graph:

        @staticmethod
        def find_cycles(graph):

            def findNewCycles(path):
                start_node = path[0]
                next_node = None
                sub = []

                # visit each edge and each node of each edge
                for edge in graph:
                    node1, node2 = edge
                    if start_node == node1:
                        next_node = node2
                    if not visited(next_node, path):
                        # neighbor node not on path yet
                        sub = [next_node]
                        sub.extend(path)
                        # explore extended path
                        findNewCycles(sub);
                    elif len(path) > 2 and next_node == path[-1]:
                        # cycle found
                        p = invert(path)
                        p = rotate_to_smallest(p)
                        if isNew(p): # and isNew(inv):
                            cycles.append(p)

            def invert(path):
                return rotate_to_smallest(path[::-1])

            #  rotate cycle path such that it begins with the smallest node
            def rotate_to_smallest(path):
                n = path.index(min(path))
                return path[n:] + path[:n]

            def isNew(path):
                return not path in cycles

            def visited(node, path):
                return node in path

            cycles = []
            for edge in graph:
                for node in edge:
                    findNewCycles([node])
            for cy in cycles:
                path = [str(node) for node in cy]
                s = ",".join(path)
            return cycles


class Convert:

    class JSON:
        @staticmethod
        def to(data):
            seq = Sequence(sequence=data['bases'])
            features = data['annotations']
            for f in features:
                F = Feature(f['name'], f['type'])
                F.color = f['color']
                seq.add_feature(f['start'], f['end'])
            seq.name = data['name']