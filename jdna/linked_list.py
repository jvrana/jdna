"""Linked list model to represent linear or circular sequences."""
from copy import copy
from functools import reduce
from typing import Any
from typing import Callable
from typing import Dict
from typing import Generator
from typing import Iterable
from typing import List
from typing import Sequence as TypingSequence
from typing import Tuple
from typing import Union


class LinkedListException(Exception):
    """Generic linked list exception."""


class LinkedListIndexError(LinkedListException, IndexError):
    """Indices were out of bounds for LinkedList."""


class Node:
    """A node in a linked list."""

    __slots__ = ["data", "__next", "__prev"]

    def __init__(self, data):
        """Node constructor. Stores a single piece of data.

        :param data: data
        :type data: any
        """
        self.data = data
        self.__next = None
        self.__prev = None

    def prev(self):
        """Return the previous node.

        :return: the previous node
        :rtype: Node
        """
        return self.__prev

    def next(self):
        """Return the next node.

        :return: the next node
        :rtype: Node
        """
        return self.__next

    def add_next(self, data):
        """Create a new node and add to next.

        :param data: any data
        :type data: any
        :return: the new node
        :rtype: Node
        """
        new_node = self.__class__(data)
        self.set_next(new_node)
        return new_node

    def add_prev(self, data):
        """Create a new node and add to previous.

        :param data: any data
        :type data: any
        :return: the new node
        :rtype: Node
        """
        new_node = self.__class__(data)
        self.set_prev(new_node)
        return new_node

    def cut_next(self):
        """Cut the next node, return the cut node.

        :return: the cut (next) node
        :rtype: Node
        """
        next_node = next(self)
        if next_node is not None:
            next_node.__prev = None
        self.__next = None
        return next_node

    def cut_prev(self):
        """Cut the previous node, return the cut node.

        :return: the cut (previous) node
        :rtype: Node
        """
        prev_node = self.prev()
        if prev_node is not None:
            prev_node.__next = None
        self.__prev = None
        return prev_node

    def _break_connections(self):
        """Break connections in this node.

        :return:
        :rtype:
        """
        self.set_next(None)
        self.set_prev(None)

    def remove(self):
        """Remove node from linked list, connecting the previous and next nodes
        together.

        :return: None
        :rtype: None
        """
        next_node = next(self)
        prev_node = self.prev()
        if next_node is not None:
            next_node.set_prev(prev_node)
        if prev_node is not None:
            prev_node.set_next(next_node)
        self._break_connections()
        return

    def swap(self):
        """Swap the previous and next nodes.

        :return: None
        :rtype: None
        """
        temp = self.__next
        self.__next = self.__prev
        self.__prev = temp

    def set_next(self, node: "Node"):
        """Set the next node.

        :param node:
        :type node:
        :return:
        :rtype:
        """
        if node is not None:
            node.__prev = self
        self.__next = node

    def set_prev(self, node: "Node"):
        """Set the previous node.

        :param node:
        :type node:
        :return:
        :rtype:
        """
        if node is not None:
            node.__next = self
        self.__prev = node

    def has_next(self) -> bool:
        return self.__next is not None

    def has_prev(self) -> bool:
        return self.__prev is not None

    def _propogate(
        self, next_method: Callable, stop: "Node" = None, stop_criteria=None
    ) -> Generator["Node", None, None]:
        visited = set()
        curr = self
        while True:
            if (
                curr is None
                or curr in visited
                or (stop_criteria and stop_criteria(curr))
            ):
                return
            yield curr
            if curr is stop:
                return
            visited.add(curr)
            curr = next_method(curr)

    def fwd(
        self, stop_node: "Node" = None, stop_criteria: Callable = None
    ) -> Generator:
        """Propogates forwards until stop node is visited or stop criteria is
        reached.

        :param stop_node:
        :type stop_node:
        :param stop_criteria:
        :type stop_criteria:
        :return:
        :rtype:
        """
        return self._propogate(
            lambda x: x.next(), stop=stop_node, stop_criteria=stop_criteria
        )

    def rev(
        self, stop_node: "Node" = None, stop_criteria: Callable = None
    ) -> Generator["Node", None, None]:
        """Propogates backwards until stop node is visited or stop criteria is
        reached.

        :param stop_node:
        :type stop_node:
        :param stop_criteria:
        :type stop_criteria:
        :return:
        :rtype:
        """
        return self._propogate(
            lambda x: x.prev(), stop=stop_node, stop_criteria=stop_criteria
        )

    def find_first(self) -> "Node":
        """Find the head node.

        :return:
        :rtype:
        """
        rev = self.rev()
        first = self
        for first in rev:
            pass
        return first

    def find_last(self) -> "Node":
        """Find the tail node.

        :return:
        :rtype:
        """
        fwd = self.fwd()
        last = self
        for last in fwd:
            pass
        return last

    def longest_match(
        self, node: "Node", next_method: Callable = None
    ) -> Tuple["Node", "Node"]:
        """Find the longest match between two linked_lists.

        :param node: the node to compare
        :type node: Node
        :param next_method: how to obtain the next node
        :type next_method: callable
        :return: list of tuples containing matching nodes
        :rtype: list
        """
        if next_method is None:

            def next_method(x):
                return lambda x: x.next()

        x1 = self
        x2 = node
        start = None
        end = None
        while x1 and x2 and x1.equivalent(x2):
            if start is None:
                start = (x1, x2)
            end = (x1, x2)
            x1 = next_method(x1)
            x2 = next_method(x2)
        return start, end

    def _complete_match(self, node: "Node", next_method: Callable) -> bool:
        """Return whether the longest match between two nodes is equivalent.

        :param node: the node to compare
        :type node: Node
        :param next_method: how to obtain the next node
        :type next_method: callable
        :return: whether the longest match between two nodes is equivalent
        :rtype: bool
        """
        length = self.longest_match(node, next_method)
        if not length:
            return False
        t1, t2 = length[-1]
        return not (next_method(t1) and next_method(t2))

    def complete_match_fwd(self, y: "Node") -> bool:
        return self._complete_match(y, lambda x: next(x))

    def complete_match_rev(self, y: "Node") -> bool:
        return self._complete_match(y, lambda x: x.prev())

    def equivalent(self, other: "Node") -> bool:
        """Evaluates whether two nodes hold the same data."""
        return self.data == other.data

    def __next__(self) -> "Node":
        return self.next()

    def copy(self) -> "Node":
        return self.__class__(self.data)

    def __copy__(self) -> "Node":
        return self.copy()

    def __deepcopy__(self, memo: Dict) -> "Node":
        raise NotImplementedError(
            "copy.deepcopy not implemented with class"
            "{}. Use copy.copy instead.".format(self.__class__.__name__)
        )

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return str(self.data)


class EmptyNode(Node):
    def __init__(self):
        pass

    # @property
    # def __next(self):
    #     return EmptyNode()
    #
    # @property
    # def __prev(self):
    #     return EmptyNode()

    def __eq__(self, other):
        return False

    def _propogate(self, *args, **kwargs):
        return
        yield

    def __str__(self):
        return ""


class LinkedListMatch:
    """A match object."""

    def __init__(
        self,
        template_bounds: Tuple[Node, Node],
        query_bounds: Tuple[Node, Node],
        template: "DoubleLinkedList",
        query: "DoubleLinkedList",
    ):
        self.span = tuple(template.indices_of(template_bounds))
        self.query_span = tuple(query.indices_of(query_bounds))
        self.query_bounds = query_bounds
        self.template_bounds = template_bounds

    @classmethod
    def batch_create(
        cls,
        template_bounds_list: List[Tuple[Node, Node]],
        query_bounds_list: List[Tuple[Node, Node]],
        template: "DoubleLinkedList",
        query: "DoubleLinkedList",
    ) -> List["LinkedListMatch"]:
        """Efficiently create several LinkedListMatches from lists of template
        starts/ends and query starts/ends.

        :param template_bounds_list: list of 2 len tuples containing starts and ends from a template
        :type template_bounds_list: template DoubleLinkedList
        :param query_bounds_list: list of 2 len tuples containing starts and ends from a query
        :type query_bounds_list: query DoubleLinkedList
        :param template: the template
        :type template: DoubleLinkedList
        :param query: the query
        :type query: DoubleLinkedList
        :return: matchese
        :rtype: list of LinkedListMatch
        """
        if len(template_bounds_list) != len(query_bounds_list):
            raise LinkedListIndexError(
                "Cannot create Matches The template bounds list must be same size as query_"
                "bounds list"
            )
        if not template_bounds_list:
            return []
        template_nodes = reduce(lambda x, y: list(x) + list(y), template_bounds_list)
        query_nodes = reduce(lambda x, y: list(x) + list(y), query_bounds_list)
        template_indices = template.indices_of(template_nodes)
        query_indices = query.indices_of(query_nodes)

        matches = []
        for i in range(0, len(template_nodes), 2):
            new_match = cls.__new__(cls)
            new_match.span = (template_indices[i], template_indices[i + 1])
            new_match.query_span = (query_indices[i], query_indices[i + 1])
            new_match.query_bounds = (query_nodes[i], query_nodes[i + 1])
            new_match.template_bounds = (template_nodes[i], template_nodes[i + 1])
            matches.append(new_match)
        return matches

    @property
    def start(self) -> Node:
        return self.template_bounds[0]

    @property
    def end(self) -> Node:
        return self.template_bounds[1]

    @property
    def query_start(self) -> Node:
        return self.query_bounds[0]

    @property
    def query_end(self) -> Node:
        return self.query_bounds[1]

    def __repr__(self) -> str:
        return "<{cls} span={span}, qspan={query_span}>".format(
            cls=self.__class__.__name__, span=self.span, query_span=self.query_span
        )

    def __str__(self) -> str:
        return self.__repr__()


class DoubleLinkedList:
    """A generic double linked list class."""

    class Direction:
        FORWARD = 1
        REVERSE = -1

    NODE_CLASS = Node

    def __init__(
        self, data: TypingSequence[Any] = None, first: Node = None, cyclic: bool = False
    ):
        """linked list construction.

        :param data: iterable data
        :type data: iterable
        :param first: first node
        :type first: Node
        """
        self._head = EmptyNode()
        if data is not None:
            self.initialize(data)
        elif first is not None:
            self._head = first
        if cyclic:
            self.circularize()

    def new_node(self, data) -> "NODE_CLASS":
        return self.NODE_CLASS(data)

    def initialize(self, sequence: TypingSequence[Any]):
        prev = None
        for i, d in enumerate(sequence):
            curr = self.new_node(d)
            if i == 0:
                self._head = curr
            if prev:
                prev.set_next(curr)
            prev = curr

    def is_empty(self) -> bool:
        return isinstance(self._head, EmptyNode)

    @property
    def head(self) -> Node:
        return self._head
        if self.is_empty():
            return self._head
        if self.cyclic:
            return self._head
        first = self._head.find_first()
        self._head = first
        return self._head

    @head.setter
    def head(self, new_head: Node) -> Node:
        self._head = new_head
        return new_head

    @property
    def tail(self) -> Node:
        return self.head.find_last()

    # TODO: This method is inefficient, but can probably be managed more manually (i.e. anytime a manipulation occurs)
    @property
    def cyclic(self) -> bool:
        if self.is_empty():
            return False
        visited = set()
        curr = self._head
        while curr:
            if curr in visited:
                return True
            visited.add(curr)
            curr = next(curr)
        return False

    @cyclic.setter
    def cyclic(self, b: bool):
        if self.cyclic and not b:
            return self.linearize()
        elif not self.cyclic and b:
            return self.circularize()

    @property
    def circular(self) -> bool:
        """Alias for cyclic."""
        return self.cyclic

    @circular.setter
    def circular(self, v: bool):
        """Alias for cyclic."""
        self.cyclic = v

    def circularize(self) -> "DoubleLinkedList":
        if not self.cyclic:
            self.tail.set_next(self.head)
        return self

    def linearize(self, i=0):
        this_i = self.get(i)
        this_i.cut_prev()
        self.head = this_i
        return self

    @property
    def nodes(self) -> List[Node]:
        n = list(self.head.fwd())
        return n

    def get(self, i: int) -> Node:
        if i is None:
            return None
        elif i < 0:
            return self.nodes[i]
        for index, n in enumerate(self):
            if index == i:
                return n
        raise LinkedListIndexError(
            "There is no node at index '{}'. There are {} nodes.".format(i, len(self))
        )

    def cut(
        self, i: Union[List[int], int, Tuple[int, ...]], cut_prev: bool = True
    ) -> List["DoubleLinkedList"]:
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
        self._check_if_in_bounds(i)
        self_copy = copy(self)
        all_nodes = self_copy.nodes
        cut_nodes = []
        for cut_loc in i:
            node = all_nodes[cut_loc]
            if cut_prev:
                c = node.cut_prev()
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
    def collect_nodes(nodes):
        """Return all visisted nodes and return an unordered set of nodes.

        :return: all visited nodes
        :rtype: set
        """
        visited = set()
        for n in nodes:
            if n not in visited:
                for tail in n.fwd(stop_criteria=lambda x: x in visited):
                    visited.add(tail)
                for head in n.rev(stop_criteria=lambda x: x in visited):
                    visited.add(head)
        return visited

    def all_nodes(self):
        return self.collect_nodes([self.head])

    @staticmethod
    def find_ends(nodes):
        """Efficiently finds the head and tails from a group of nodes."""
        visited = set()
        heads = []
        tails = []
        for n in nodes:
            if n not in visited:
                head = n
                tail = n
                visited_tails = set()
                visited_heads = set()
                for tail in n.fwd(stop_criteria=lambda x: x in visited):
                    visited_tails.add(tail)
                for head in n.rev(stop_criteria=lambda x: x in visited):
                    visited_heads.add(head)
                visited = visited.union(visited_heads)
                visited = visited.union(visited_tails)
                if head not in heads and tails not in tails:
                    if not head.has_prev() is None and not tail.has_next():
                        heads.append(head)
                        tails.append(tail)
        return zip(heads, tails)

    @classmethod
    def segments(cls, nodes: Iterable[Node]) -> List["DoubleLinkedList"]:
        return [cls(first=h) for h, _ in cls.find_ends(nodes)]

    def insert(
        self, node_list: List[Node], i: int, copy_insertion: bool = True
    ) -> "DoubleLinkedList":
        if i == len(self):
            pass
        else:
            self._check_if_in_bounds(i)
        if node_list.cyclic:
            raise TypeError("Cannot insert a cyclic sequence")
        if copy_insertion:
            node_list = copy(node_list)
        # TODO: This copies the insertion sequence, you want that?
        if i == len(self):
            loc2 = None
            loc1 = self.get(i - 1)
        else:
            loc2 = self.get(i)
            loc1 = loc2.prev()
        first = node_list.nodes[0]
        last = node_list.nodes[-1]
        first.set_prev(loc1)
        last.set_next(loc2)
        if (
            i == 0
        ):  # Special case in which user inserts sequence in front of their sequence; they probably intend to re-index it
            self.head = first
        return self

    def remove(self, i: int) -> "DoubleLinkedList":
        self._check_if_in_bounds(i)
        to_be_removed = self.get(i)
        new_first = self.head
        if i == 0:
            new_first = next(new_first)
        to_be_removed.remove()
        self.head = new_first
        return self

    def reindex(self, i: int) -> "DoubleLinkedList":
        self._check_if_in_bounds(i)
        if not self.cyclic:
            raise TypeError(
                "Cannot re-index a linear {}".format(self.__class__.__name__)
            )
        self.head = self.get(i)
        return self

    def _check_if_in_bounds(self, num: int):
        if isinstance(num, int):
            num = [num]
        for n in num:
            mn = 0
            mx = len(self) - 1
            if n < 0 or n > mx:
                raise LinkedListIndexError(
                    "Index {} out of acceptable bounds ({}, {})".format(n, mn, mx)
                )

    # TODO: implement yield in find_iter, search_all should call this
    # TODO: query should be any interable
    # TODO: [:1] and [-10:] style cuts should be available
    # TODO: documentation for methods
    # TODO: move DoubleLinkedList to its own thing?
    # TODO: element insertion
    # TODO: search should return a 'cut' of the sequence
    def search_all(self, query: "DoubleLinkedList") -> List[Tuple[int, Node]]:
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

    # def longest_match(self, other):
    #     n1 = self.head
    #     n2 = other.head
    #     start = None
    #     while n1 and n2 and n1.equivalent(n2):
    #         if start is None:
    #             start = (n1, n2)
    #         n1 = next(n1)
    #         n2 = next(n2)
    #
    #         curr = next(curr)
    #     for n in self:
    #         if other_node.equivalent(n):
    #             stop = n
    #             if start is None:
    #                 start = stop
    #         other_node = next(other_node)

    # def longest_match(self, other):
    #     matched = []
    #     i = 0
    #     j = i
    #     for x1, x2 in zip(self, other):
    #         if x1.equivalent(x2):
    #             j += 1
    #             matched.append((x1, x2))
    #     if j > len(self):
    #         if not self.cyclic:
    #             raise Exception("Template was expected to be cyclic but was not")
    #         j -= len(self)
    #     return LinkedListMatch(self.head, )
    #     if len(matched) == qlen:
    #         if j > len(self):
    #             if self.cyclic:
    #                 j -= len(self)
    #             else:
    #                 raise Exception("Template was expected to be cyclic but was not")
    #         yield LinkedListMatch(curr_node, x1, span=(i, j - 1))

    @classmethod
    def match(
        cls,
        n1: Node,
        n2: Node,
        query_direction: int = Direction.FORWARD,
        template_direction: int = Direction.FORWARD,
        protocol: Callable = None,
    ):
        """"""
        if protocol is None:

            def protocol(x, y):
                return x.equivalent(y)

        iterators = {
            cls.Direction.FORWARD: lambda x: x.fwd(),
            cls.Direction.REVERSE: lambda x: x.rev(),
        }
        for x1, x2 in zip(
            iterators[template_direction](n1), iterators[query_direction](n2)
        ):
            if protocol(x1, x2):
                yield (x1, x2)
            else:
                return

    def find_iter(
        self,
        query: "DoubleLinkedList",
        min_query_length: int = None,
        direction: int = Direction.FORWARD,
        protocol: Callable = None,
        depth: int = None,
    ):
        """Iteratively finds positions that match the query.

        :param query: query list to find
        :type query: DoubleLinkedList
        :param min_query_length: the minimum number of matches to return. If None, find_iter will only return complete
                                matches
        :type min_query_length: inst
        :param direction: If Direction.FORWARD (+1), find iter will search from the query head and search forward,
                        potentially leaving a 'tail' overhang on the query. If
                        Direction.REVERSE (-1), find iter will search from the query tail and search reverse, potentially
                        leaving a 'head' overhang on the query. If a tuple, the template_direction and query_direction are
                        set respectively.
        :type direction: int or tuple
        :param protocol: the callable taking two parameters (as Node) to compare during find. If None, defaults
                        to 'equivalent'
        :type protocol: callable
        :return: list of LinkedListMatches
        :rtype: list
        """
        if isinstance(direction, tuple):
            template_direction, query_direction = direction
        else:
            template_direction, query_direction = direction, direction

        curr_node = self.head
        visited = set()
        if self.Direction.REVERSE == query_direction:
            query_start = query.tail
        else:
            query_start = query.head
        qlen = len(query)
        if min_query_length is None:
            min_query_length = qlen

        template_nodes = []
        query_nodes = []
        index = 0
        while curr_node and curr_node not in visited:
            index += 1
            if depth is not None and index > depth:
                break
            visited.add(curr_node)
            matches = list(
                self.match(
                    curr_node,
                    query_start,
                    query_direction=query_direction,
                    template_direction=template_direction,
                    protocol=protocol,
                )
            )
            if self.Direction.REVERSE == template_direction:
                matches = matches[::-1]
            if len(matches) >= min_query_length:
                query_nodes.append([matches[0][1], matches[-1][1]])
                template_nodes.append([matches[0][0], matches[-1][0]])
            curr_node = next(curr_node)
        return LinkedListMatch.batch_create(template_nodes, query_nodes, self, query)

    def reverse(self) -> "DoubleLinkedList":
        if self.is_empty():
            return self
        nodes = self.nodes
        for s in nodes:
            s.swap()
        self.head = nodes[-1]
        # if self.cyclic:
        #     self.reindex(1)
        return self

    def fuse_in_place(self, seq: "DoubleLinkedList") -> "DoubleLinkedList":
        if seq.is_empty():
            return self
        elif self.is_empty():
            return seq.copy()
        f = self.head.find_last()
        length = seq.head
        f.set_next(length)
        return self

    def fuse(self, other: "DoubleLinkedList") -> "DoubleLinkedList":
        return copy(self).fuse_in_place(copy(other))

    def copy(self) -> "DoubleLinkedList":
        return self.__copy__()

    @classmethod
    def empty(cls) -> "DoubleLinkedList":
        return cls()

    @staticmethod
    def empty_iterator() -> Generator:
        return
        yield

    def node_range(self, start: Node, end: Node) -> Generator[Node, None, None]:
        """Iterate between 'start' to 'end' nodes (inclusive)"""
        for n in start.fwd(stop_node=end):
            yield n

    @classmethod
    def new_slice(cls, start: Node, end: Node) -> "DoubleLinkedList":
        """Return a copy of the sequence between 'start' and 'end' nodes)"""
        nodes = []
        if start is not None:
            for n in start.fwd(stop_node=end):
                nodes.append(n)
                # new_node = copy(n)
                # new_nodes.append(new_node)
                # new_node.set_prev(prev)
                # prev = new_node
        else:
            if end is None:
                return cls.empty()
            for n in end.rev(stop_node=start):
                nodes.append(n)
            nodes = nodes[::-1]
        if end is not None and nodes[-1] is not end:
            return cls()
        if start is not None and nodes[0] is not start:
            return cls()

        next = None
        for n in nodes[::-1]:
            n = copy(n)
            n.set_next(next)
            next = n
        n.set_prev(None)

        return cls(first=n)

    def copy_slice(self, start: Node, end: Node) -> "DoubleLinkedList":
        """Return a copy of the sequence between 'start' and 'end' nodes.

        If start is None, return the slice copy from head to end. If end
        is None, return copy from start to tail. If both start and end
        are None return None.
        """
        if start is None:
            start = self.head
        return self.new_slice(start, end)

    def range(self, i: int, j: int) -> Generator[Node, None, None]:
        """Returns an iterator from node at 'i' to node at 'j-1'."""
        if i == j:
            return self.empty_iterator()
        if j is not None:
            if j == 0 and not self.cyclic:
                return self.empty_iterator()
            j = j - 1
        try:
            return self.inclusive_range(i, j)
        except LinkedListIndexError:
            return self.empty_iterator()

    def inclusive_range(self, i: int, j: int) -> Generator[Node, None, None]:
        """Return generator for inclusive nodes between index i and j.

        If i is None, assume i is the head node.
        """
        start = self.get(i)
        if start is None:
            start = self.head
        try:
            stop = self.get(j)
        except LinkedListIndexError:
            stop = None
        stop_hit = False
        for n in start.fwd(stop_node=stop):
            yield n
            if n is stop:
                stop_hit = True
        if stop is not None and not stop_hit:
            raise LinkedListIndexError(
                "Inclusive indices {} out of bounds".format((i, j))
            )

    def index_of(self, node: Node) -> int:
        for i, n in enumerate(self):
            if n is node:
                return i

    def indices_of(self, nodes: Iterable[Node]) -> List[int]:
        index_dict = {}
        for i, n in enumerate(self):
            for node in nodes:
                if n is node:
                    index_dict[n] = i
        return [index_dict.get(n, None) for n in nodes]

    def data(self) -> List[Any]:
        return [n.data for n in self]

    def compare(self, other: "DoubleLinkedList") -> bool:
        """Compares two linked lists. If both are cyclic, will attempt to
        reindex.

        :param other: other linked list
        :type other: DoubleLinkedList
        :return: whether sequence data is equivalent
        :rtype: bool
        """
        if self.cyclic and other.cyclic:
            anchors = []
            if len(other) > 100:
                for a in self.find_iter(other[:20]):
                    anchors.append(a.span[0])
            else:
                anchors = range(len(other))
            for a in anchors:
                temp = self.copy().reindex(a)
                if temp == other:
                    return True
            return False
        else:
            return self == other

    def left_trim(self, i: int) -> "DoubleLinkedList":
        if self.cyclic:
            raise IndexError("Cannot chop a cyclic sequence.")
        return self.cut(i)[-1]

    def right_trim(self, i: int) -> "DoubleLinkedList":
        if self.cyclic:
            raise IndexError("Cannot chop a cyclic sequence.")
        return self.cut(i, cut_prev=False)[0]

    def __eq__(self, other: "DoubleLinkedList") -> bool:
        if self.cyclic != other.cyclic:
            return False
        data1 = [n.data for n in self]
        data2 = [n.data for n in other]
        return data1 == data2

    def __add__(self, other: "DoubleLinkedList") -> "DoubleLinkedList":
        return self.fuse(other)

    def __getitem__(self, key: Union[slice, int]) -> Union["DoubleLinkedList", Node]:
        if isinstance(key, slice):
            if key.step and key.step > 1:
                raise LinkedListIndexError(
                    "Step > 1 is not supported for sliced object of '{}'".format(
                        self.__class__.__name__
                    )
                )

            if key.start is not None and key.stop is not None:
                if key.start == key.stop:
                    return self.empty()

            new_list = self.__copy__()
            if key.start is None and key.stop is None:
                if key.step == -1:
                    return new_list.reverse()
                else:
                    return new_list

            try:
                new_nodes = list(new_list.range(key.start, key.stop))
            except LinkedListIndexError:
                return self.empty()
            if not new_nodes:
                return self.empty()
            start = new_nodes[0]
            end = new_nodes[-1]

            start.cut_prev()
            end.cut_next()
            return self.__class__(first=start)
        return self.get(key)

    def __contains__(self, item: Node) -> bool:
        return item in self.all_nodes()

    def __copy__(self) -> "DoubleLinkedList":
        head = None
        prev = None
        for n in self.nodes:
            new_node = n.copy()
            if head is None:
                head = new_node
            new_node.set_prev(prev)
            prev = new_node
        copied = type(self)(first=head)
        if self.cyclic:
            copied.cyclic = self.cyclic
        return copied

    def __deepcopy__(self, memo: Dict) -> "DoubleLinkedList":
        raise NotImplementedError(
            "copy.deepcopy not implemented with class"
            "{}. Use copy.copy instead.".format(self.__class__.__name__)
        )

    def __reversed__(self) -> "DoubleLinkedList":
        copied = self.copy()
        copied.reverse()
        return copied

    def __len__(self) -> int:
        length = 0
        for _ in self:
            length += 1
        return length

    def __iter__(self) -> Generator[Node, None, None]:
        return self.head.fwd()

    def __repr__(self) -> str:
        return "<{cls} data='{data}'>".format(
            cls=self.__class__.__name__, data=str(self)
        )
        return str(self)

    def __str__(self) -> str:
        return "".join(str(x) for x in self)
