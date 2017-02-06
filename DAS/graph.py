graph = {
    1: [2, 3, 4],
    2: [5, 6],
    3: [10],
    4: [7, 8],
    5: [9, 10],
    7: [11, 12],
    11: [13],
    8: ['end']
}


def bfs(graph, queue, paths):
    while queue:
        path = queue.pop(0)
        if path[-1] in graph:
            print '\t', path
            for node in graph[path[-1]]:
                new_path = path[:] + [node]
                queue.append(new_path)
                bfs(graph, queue, paths)
        else:
            paths.append(path[:])
    return paths

def dfs_iter(graph, root):
    paths = []
    stack = [[root]]

    while stack:

        path = stack.pop()
        n = path[-1]
        if n in graph:
            for node in graph[n]:
                new_path = path[:] + [node]
                stack.append(new_path)
        else:
            paths.append(path)
    return paths


print dfs_iter(graph, 1)
