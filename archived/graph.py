"""
graph utilities for reactions
"""

import itertools

class Graph:

    @classmethod
    def find_cycles(cls, graph, min_path_length=2):

        unique_cycles = []
        for g in graph:
            cycles = cls.find_cyclic_paths(graph, [g], min_path_length=min_path_length)
            for cy in cycles:
                cy = cls.rotate_to_smallest(cy)
                if cy not in unique_cycles:
                    unique_cycles.append(cy)
        return unique_cycles

    @classmethod
    def find_linear(cls, graph, min_path_length=2):
        unique_paths = []
        for g in graph:
            c = True
            for p in unique_paths:
                if g in p:
                    c = False
                    break
            if c:
                unique_paths += cls.find_linear_paths(graph, [g], min_path_length=min_path_length)
        return unique_paths

    @classmethod
    def find_linear_paths(cls, graph, path, min_path_length=2):
        if path[-1] not in graph and len(path) >= min_path_length:
            return [path]
        paths = []
        for node in graph[path[-1]]:
            if node not in path:
                newpaths = cls.find_linear_paths(graph, path[:] + [node],
                                                             min_path_length=min_path_length)
                paths += newpaths
        return paths

    @classmethod
    def find_cyclic_paths(cls, graph, path, min_path_length=2):
        if len(path) >= min_path_length and path[-1] == path[0]:
            return [path[:-1]]
        paths = []
        if path[-1] not in graph:
            return []
        for node in graph[path[-1]]:
            if node not in path[1:]:
                newpaths = cls.find_cyclic_paths(graph, path[:] + [node],
                                                             min_path_length=min_path_length)
                paths += newpaths
        return paths


    @staticmethod
    def rotate_to_smallest(path):
        n = path.index(min(path))
        return path[n:] + path[:n]


    @staticmethod
    def group_ranges(lst):
        pos = (j - i for i, j in enumerate(lst))
        t = 0
        for i, els in itertools.groupby(pos):
            l = len(list(els))
            el = lst[t]
            t += l
            yield list(range(el, el + l))


    @staticmethod
    def lcs(S, T):
        m = len(S)
        n = len(T)
        counter = [[0] * (n + 1) for x in range(m + 1)]
        longest = 0
        lcs_set = set()
        for i in range(m):
            for j in range(n):
                if S[i] == T[j]:
                    c = counter[i][j] + 1
                    counter[i + 1][j + 1] = c
                    if c > longest:
                        lcs_set = set()
                        longest = c
                        lcs_set.add(S[i - c + 1:i + 1])
                    elif c == longest:
                        lcs_set.add(S[i - c + 1:i + 1])

        return lcs_set