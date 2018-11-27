"""
Simulate molecular reactions
"""

# import itertools
# import warnings
# from collections import defaultdict
# from copy import copy
# from jdna.sequence import Sequence
# from jdna.graph import Graph

from jdna.sequence import Sequence, SequenceFlags
from collections import namedtuple
import networkx as nx
import itertools

class ReactionException(Exception):
    """Generic reaction exception"""


class PCRException(ReactionException):
    """Exception with pcr"""

BindingEvent = namedtuple("BindingEvent", ["template", "primer", "position"])

class Reaction(object):
    MIN_BASES = 13

    @classmethod
    def pcr(cls, template, p1, p2, min_bases=MIN_BASES):
        bindings = list(template.anneal(p1, min_bases=min_bases))
        bindings += template.anneal(p2, min_bases=min_bases)

        forward_bindings = []
        reverse_bindings = []
        for bind in bindings:
            if bind.direction == SequenceFlags.FORWARD:
                forward_bindings.append(bind)
            elif bind.direction == SequenceFlags.REVERSE:
                reverse_bindings.append(bind)
            else:
                raise Exception("Direction {} not recognized".format(bind.direction))

        if not forward_bindings or not reverse_bindings:
            raise PCRException("Some primers did not bind. Number of forward bindings: {}. Number of rev bindings: {}".format(
                len(forward_bindings), len(reverse_bindings)
            ))
        products = []
        for fwd, rev in itertools.product(forward_bindings, reverse_bindings):
            overhang1 = fwd.five_prime_overhang
            overhang2 = rev.five_prime_overhang
            print(fwd.span)
            print(rev.span)
            print(template.index_of(fwd.start))
            print(template.index_of(rev.end))
            amplified = template.new_slice(fwd.start, rev.end)
            print(overhang1)
            print(overhang2)
            products.append(overhang1 + amplified + overhang2)
        return products

    @classmethod
    def anneal_sequences(cls, sequences, min_bases=Sequence.DEFAULTS.MIN_ANNEAL_BASES):
        # pairs = itertools.product(sequences, sequences + [s.copy().reverse_complement() for s in sequences[1:]])
        pairs = itertools.product(sequences, repeat=2)
        bindings = []
        for s1, s2 in pairs:
            for binding in s1.dsanneal(s2, min_bases=min_bases):
                if not (binding.span[0] == 0 and binding.span[-1] == len(s1)-1):
                    bindings.append(BindingEvent(s1, s2, binding))
        return bindings

    @classmethod
    def interaction_graph(cls, sequences, min_bases=Sequence.DEFAULTS.MIN_ANNEAL_BASES, bind_reverse_complement=False):
        G = nx.DiGraph()
        for s in sequences:
            G.add_node(s.global_id, sequence=s)
        if bind_reverse_complement:
            sequences = [s.copy() for s in sequences] + [s.copy().reverse_complement() for s in sequences]
        seq_pairs = itertools.product(sequences, repeat=2)
        for s1, s2 in seq_pairs:
            for binding in s1.anneal_forward(s2):
                if not (binding.span[0] == 0 and binding.span[-1] == len(s1) - 1):
                    G.add_edge(s2.global_id, s1.global_id, template=s1, primer=s2, binding=binding)
        return G

    @classmethod
    def linear_paths(cls, G):
        paths = []
        subgraph = G.subgraph(G.nodes).copy()
        for _ in nx.simple_cycles(G):
            return []
        while subgraph.nodes:
            path = nx.dag_longest_path(subgraph)
            if len(path) < 2:
                break
            subgraph.remove_nodes_from(path)
            paths.append(path)
        return paths

    @classmethod
    def cyclic_paths(cls, G):
        paths = []
        for path in nx.simple_cycles(G.copy()):
            paths.append(path)
        return paths

    @classmethod
    def path_to_edge_data(cls, G, path, cyclic):
        if cyclic:
            path_pairs = zip(path, path + [path[0]])
        else:
            path_pairs = zip(path[:-1], path[1:])
        edges = []
        for n1, n2 in path_pairs:
            edges.append(G.edges[n1, n2])
        return edges

    @classmethod
    def linear_assemblies(cls, sequences, min_bases=Sequence.DEFAULTS.MIN_ANNEAL_BASES):
        G = cls.interaction_graph(sequences, bind_reverse_complement=True, min_bases=min_bases)
        linear_paths = cls.linear_paths(G)
        if not linear_paths:
            return []

        for path in linear_paths:
            prev_query_end = None
            node_pairs = []
            overhangs = []
            for i, data in enumerate(cls.path_to_edge_data(G, path, False)):
                binding = data['binding']
                node_pairs.append((prev_query_end, binding.query_start))
                prev_query_end = binding.query_end
                overhangs.append(binding.template_anneal)
            node_pairs.append((prev_query_end, None))
            overhangs.append(Sequence())

            amplified_sequences = [Sequence.new_slice(*pair) for pair in node_pairs]

            import functools
            product = functools.reduce(lambda x, y: x + y, zip(amplified_sequences, overhangs))
            print(len(product))
