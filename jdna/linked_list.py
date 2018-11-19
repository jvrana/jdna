from copy import copy


class Node(object):
    """
    A node in a node_ list
    """

    def __init__(self, data):
        self.data = data
        self.__next = None
        self.__prev = None

    def __next__(self):
        return self.__next

    @property
    def prev(self):
        """
        Return the previous node

        :return: the previous node
        :rtype: Node
        """
        return self.__prev

    @property
    def next(self):
        """
        Return the next node

        :return: the next node
        :rtype: Node
        """
        return self.__next

    def add_next(self, data):
        """
        Create a new node and add to next

        :param data: any data
        :type data: any
        :return: the new node
        :rtype: Node
        """
        new_node = Node(data)
        self.set_next(new_node)
        return new_node

    def add_prev(self, data):
        """
        Create a new node and add to previous

        :param data: any data
        :type data: any
        :return: the new node
        :rtype: Node
        """
        new_node = Node(data)
        self.set_prev(new_node)
        return new_node

    def cut_next(self):
        """
        Cut the next node, return the cut node.

        :return: the cut (next) node
        :rtype: Node
        """
        next_node = next(self)
        if next_node is not None:
            next_node.__assign_prev(None)
        self.__assign_next(None)
        return next_node

    def cut_prev(self):
        """
        Cut the previous node, return the cut node.

        :return: the cut (previous) node
        :rtype: Node
        """
        prev_node = self.prev
        if prev_node is not None:
            prev_node.__assign_next(None)
        self.__assign_prev(None)
        return prev_node

    def _break_connections(self):
        """
        Break connections in this node.
        :return:
        :rtype:
        """
        self.set_next(None)
        self.set_prev(None)

    def remove(self):
        """
        Remove node from linked list, connecting the previous and next nodes together.

        :return: None
        :rtype: None
        """
        next_node = next(self)
        prev_node = self.prev
        if next_node is not None:
            next_node.set_prev(prev_node)
        if prev_node is not None:
            prev_node.set_next(next_node)
        self._break_connections()
        return

    def swap(self):
        """
        Swap the previous and next nodes.

        :return: None
        :rtype: None
        """
        temp = self.__next
        self.__next = self.__prev
        self.__prev = temp

    def __assign_next(self, node):
        self.__next = node

    def __assign_prev(self, node):
        self.__prev = node

    def set_next(self, node):
        """
        Set the next node

        :param node:
        :type node:
        :return:
        :rtype:
        """
        if node is not None:
            node.__assign_prev(self)
        self.__assign_next(node)

    def set_prev(self, node):
        """
        Set the previous node

        :param node:
        :type node:
        :return:
        :rtype:
        """
        if node is not None:
            node.__assign_next(self)
        self.__assign_prev(node)

    def _propogate(self, next_method, stop=None, stop_criteria=None):
        visited = set()
        curr = self
        while True:
            if curr is None or curr in visited or (stop_criteria and stop_criteria(curr)):
                return
            yield curr
            visited.add(curr)
            if curr is stop:
                return
            curr = next_method(curr)

    def fwd(self, stop_node=None, stop_criteria=None):
        """
        Propogates forwards until stop node is visited or stop criteria is reached.

        :param stop_node:
        :type stop_node:
        :param stop_criteria:
        :type stop_criteria:
        :return:
        :rtype:
        """
        return self._propogate(
            lambda x: next(x),
            stop=stop_node,
            stop_criteria=stop_criteria)

    def rev(self, stop_node=None, stop_criteria=None):
        """
        Propogates backwards until stop node is visited or stop criteria is reached.

        :param stop_node:
        :type stop_node:
        :param stop_criteria:
        :type stop_criteria:
        :return:
        :rtype:
        """
        return self._propogate(
            lambda x: x.prev,
            stop=stop_node,
            stop_criteria=stop_criteria)

    def find_first(self):
        """
        Find the head node

        :return:
        :rtype:
        """
        rev = self.rev()
        first = self
        for n in rev:
            first = n
        return first

    def find_last(self):
        """
        Find the tail node

        :return:
        :rtype:
        """
        fwd = self.fwd()
        while True:
            try:
                n = next(fwd)
            except StopIteration:
                return n

    def _longest_match(self, node, next_method):
        """
        Find the longest match between two linked_lists

        :param node: the node to compare
        :type node: Node
        :param next_method: how to obtain the next node
        :type next_method: callable
        :return: list of tuples containing matching nodes
        :rtype: list
        """
        x1 = self
        x2 = node
        longest_match = []
        while x1 and x2:
            if x1.equivalent(x2):
                longest_match.append((x1, x2))
                x1 = next_method(x1)
                x2 = next_method(x2)
            else:
                break
        return longest_match

    def _complete_match(self, node, next_method):
        """
        Return whether the longest match between two nodes is equivalent.

        :param node: the node to compare
        :type node: Node
        :param next_method: how to obtain the next node
        :type next_method: callable
        :return: whether the longest match between two nodes is equivalent
        :rtype: bool
        """
        l = self._longest_match(node, next_method)
        if not l:
            return False
        t1, t2 = l[-1]
        return not (next_method(t1) and next_method(t2))

    def complete_match_fwd(self, y):
        return self._complete_match(y, lambda x: next(x))

    def complete_match_rev(self, y):
        return self._complete_match(y, lambda x: x.prev)

    def equivalent(self, other):
        """Evaluates whether two nodes hold the same data"""
        return self.data == other.data

    def __copy__(self):
        copied = type(self)(self.data)
        return copied

    def __deepcopy__(self, memo):
        raise NotImplementedError("copy.deepcopy not implemented with class" \
                                  "{}. Use copy.copy instead.".format(self.__class__.__name__))

    def __repr__(self):
        return str(self)

    def __str__(self):
        return str(self.data)


class DoubleLinkedList(object):

    NODE_CLASS = Node

    def __init__(self, first=None, sequence=None):
        self._head = None
        self.not_head = None
        if sequence is not None:
            self.initialize(sequence)
            print("HEAD: {}".format(self._head))
        elif first is not None:
            self._head = first

    # @classmethod
    # def initialize_by_node(cls):

    # def find(self, data):
    #     t = self.head
    #     visited = set()
    #     while t and t not in visited:
    #         next_method = lambda x: next(x)

    def new_node(self, data):
        return self.NODE_CLASS(data)

    def initialize(self, sequence):
        prev = None
        for i, d in enumerate(sequence):
            curr = self.new_node(d)
            if i == 0:
                self._head = curr
                self.not_head = curr
            if prev:
                prev.set_next(curr)
            prev = curr

    @property
    def head(self):
        if self.cyclic:
            return self._head
        assert self._head
        first = self._head.find_first()
        assert first, "First cannot be None"
        self._head = first
        return self._head
    #
    # @head.setter
    # def head(self, node):
    #     self._head = node

    @property
    def tail(self):
        return self.head.find_last()

    @property
    def cyclic(self):
        visited = set()
        curr = self._head
        while curr:
            if curr in visited:
                return True
            visited.add(curr)
            curr = next(curr)
        return False

    @cyclic.setter
    def cyclic(self, b):
        if self.cyclic and not b:
            return self.linearize()
        elif not self.cyclic and b:
            return self.circularize()

    def circularize(self):
        if not self.cyclic:
            return self.tail.set_next(self.head)

    def linearize(self, i=0):
        this_i = self.nodes[i]
        this_i.cut_prev
        return this_i

    @property
    def nodes(self):
        return list(self.head.fwd())

    def get(self, i):
        index = 0
        for n in self:
            if index == i:
                return n
            index += 1
        raise IndexError("There is no node at index '{}'. There are only {} nodes.".format(i, index))

    def cut(self, i, cut_prev=True):
        if isinstance(i, tuple):
            i = list(i)
        if isinstance(i, int):
            i = [i]
        # Special case in which i == len
        i = list(set(i))
        if len(self) in i and cut_prev:
            i.remove(len(self))
            if self.cyclic:
                i.append(0)
        i = list(set(i))
        i.sort()
        self._inbounds(i)
        self_copy = copy(self)
        all_nodes = self_copy.nodes
        cut_nodes = []
        for cut_loc in i:
            node = all_nodes[cut_loc]
            c = None
            if cut_prev:
                c = node.cut_prev
                if c is not None:
                    cut_nodes.append(c)
                cut_nodes.append(node)
            else:
                cut_nodes.append(node)
                c = node.cut_next()
                if c is not None:
                    cut_nodes.append(c)
        return self.segments(cut_nodes)

    @staticmethod
    def all_nodes(nodes):
        """
        Return all visisted nodes and return an unordered set of nodes.

        :return: all visited nodes
        :rtype: set
        """
        visited = set()
        for n in nodes:
            if n not in visited:
                for tail in n.fwd(stop_criteria=lambda x: x not in visited):
                    visited.add(tail)
                for head in n.rev(stop_criteria=lambda x: x not in visited):
                    visited.add(head)
        return visited

    @staticmethod
    def find_ends(nodes):
        """Efficiently finds the head and tails from a group of nodes."""
        visited = set()
        pairs = set()
        for n in nodes:
            if n not in visited:
                for tail in n.fwd(stop_criteria=lambda x: x not in visited):
                    visited.add(tail)
                for head in n.rev(stop_criteria=lambda x: x not in visited):
                    visited.add(head)
                pairs.add([head, tail])

    @classmethod
    def segments(cls, nodes):
        return [cls(first=h) for h, _ in cls.find_ends(nodes)]

    def insert(self, node_list, i, copy_insertion=True):
        if i == len(self.nodes):
            pass
        else:
            self._inbounds(i)
        if node_list.cyclic:
            raise TypeError("Cannot insert a cyclic sequence")
        if copy_insertion:
            node_list = copy(node_list)
        # TODO: This copies the insertion sequence, you want that?
        if i == len(self.nodes):
            loc2 = None
            loc1 = self.nodes[i - 1]
        else:
            loc2 = self.nodes[i]
            loc1 = loc2.prev
        first = node_list.nodes[0]
        last = node_list.nodes[-1]
        first.set_prev(loc1)
        last.set_next(loc2)
        if i == 0:  # Special case in which user inserts sequence in front of their sequence; they probably intend to re-index it
            self._head = first
        return self

    def remove(self, i):
        self._inbounds(i)
        to_be_removed = self.nodes[i]
        new_first = self.head
        if i == 0:
            new_first = next(new_first)
        to_be_removed.remove()
        self._head = new_first
        return

    def reindex(self, i):
        self._inbounds(i)
        if not self.cyclic:
            raise TypeError("Cannot re-index a linear node_ set")
        self._head = self.nodes[i]

    def _inbounds(self, num):
        if isinstance(num, int):
            num = [num]
        for n in num:
            mn = 0
            mx = len(self.nodes) - 1
            if n < 0 or n > mx:
                raise IndexError("Index {} out of acceptable bounds ({}, {})".format(n, mn, mx))

    # TODO: implement yield in find_iter, search_all should call this
    # TODO: query should be any interable
    # TODO: [:1] and [-10:] style cuts should be available
    # TODO: documentation for methods
    # TODO: move DoubleLinkedList to its own thing?
    # TODO: element insertion
    # TODO: search should return a 'cut' of the sequence
    def search_all(self, query):
        curr_node = self.head
        q_node = query.head
        i = 0
        found = []
        visited = set()
        while curr_node and curr_node not in visited:
            visited.add(curr_node)
            if curr_node.complete_match_fwd(q_node):
                found.append((i, curr_node))
            curr_node = next(curr_node)
            i += 1
        return found

    def reverse(self):
        for s in self.nodes:
            s.swap()
        if self.cyclic:
            self.reindex(1)
        return self

    def __getitem__(self, key):
        if isinstance(key, slice):
            new_list = self.__copy__()
            start = new_list.nodes[key.start]
            end = new_list.nodes[key.stop - 1]
            start.cut_prev
            end.cut_next()
            return self.__class__(first=start)
        return self.nodes[key].data

    def __copy__(self):
        copied = type(self)(sequence='X')
        copied.__dict__.update(self.__dict__)
        copied.initialize(str(self))
        if self.cyclic:
            copied.circularize()
        return copied

    def __deepcopy__(self, memo):
        raise NotImplementedError("copy.deepcopy not implemented with class" \
                                  "{}. Use copy.copy instead.".format(self.__class__.__name__))

    def __reversed__(self):
        for s in self.nodes:
            s.swap()
        return self

    def __len__(self):
        return len(self.nodes)

    def __iter__(self):
        current = self._head
        visited = set()
        while current is not None:
            if current in visited:
                raise StopIteration
            yield current
            visited.add(current)
            current = next(current)

    def __repr__(self):
        return "<{cls} data='{data}'>".format(
            self.__class__.__name__,
            str(self)
        )
        return str(self)

    def __str__(self):
        return ''.join(str(x) for x in self.nodes)
