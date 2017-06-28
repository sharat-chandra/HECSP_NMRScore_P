def find_all_paths(graph, start, end, path=[]):
        path = path + [start]
        if start == end:
            return [path]
        if start not in graph:
            return []
        paths = []
        for node in graph[start]:
            if node not in path:
                paths += find_all_paths(graph, node, end, path)
        return paths

#print(find_all_paths(graph, 'A', 'B'))
