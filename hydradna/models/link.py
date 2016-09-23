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
        if not l:
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