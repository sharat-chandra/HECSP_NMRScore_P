#graph = {'A': set(['B', 'C']),
#         'B': set(['A', 'D', 'E']),
#         'C': set(['A', 'F']),
#         'D': set(['B']),
#         'E': set(['B', 'F']),
#         'F': set(['C', 'E']),
#         'G': set(['H']),
#         'H': set(['G'])        }
#def dfs(graph, start):
#    visited, stack = set(), [start]
#    while stack:
#        vertex = stack.pop()
#        if vertex not in visited:
#            visited.add(vertex)
#            stack.extend(graph[vertex] - visited)
#    return visited
#def dfs_paths(graph, start, goal):
#    stack = [(start, [start])]
#    while stack:
#        (vertex, path) = stack.pop()
#        for next in graph[vertex] - set(path):
#            if next == goal:
#                yield path + [next]
#            else:
#                stack.append((next, path + [next]))
#
#result=list(dfs_paths(graph, 'A', 'A'))
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
