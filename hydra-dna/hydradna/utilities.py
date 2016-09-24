from itertools import groupby


class Utilities:


    class Graph:

        @staticmethod
        def find_cycles(graph, min_path_length=1):
            def find_new_cycles(path):
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
                        find_new_cycles(sub)
                    elif len(path) >= min_path_length and next_node == path[-1]:
                        # cycle found
                        p = invert(path)
                        p = rotate_to_smallest(p)
                        if is_new(p): # and isNew(inv):
                            cycles.append(p)

            def invert(pth):
                return rotate_to_smallest(pth[::-1])

            #  rotate cycle path such that it begins with the smallest node
            def rotate_to_smallest(path):
                n = path.index(min(path))
                return path[n:] + path[:n]

            def is_new(path):
                return path not in cycles

            def visited(node, path):
                return node in path

            cycles = []
            for e in graph:
                for n in e:
                    find_new_cycles([n])
            return cycles

    @staticmethod
    def group_ranges(lst):
        pos = (j - i for i, j in enumerate(lst))
        t = 0
        for i, els in groupby(pos):
            l = len(list(els))
            el = lst[t]
            t += l
            yield range(el, el + l)