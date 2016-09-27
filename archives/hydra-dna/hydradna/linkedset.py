from hydradna.link import Link

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
        if len(self) in i and cut_prev:
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
        if i == 0:  # Special case in which user inserts sequence in front of their sequence; they probably intend to re-index it
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

    def reverse(self):
        for s in self.get():
            s.swap()
        if self.is_cyclic():
            self.reindex(1)
        return self

    def copy(self):
        return deepcopy(self)

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
