
def get_connected_group(node,graph, already_seen):
        result = []
        nodes = set([node])
        while nodes:
            node = nodes.pop()
            already_seen.add(node)
            nodes = nodes or graph[node] - already_seen
            result.append(node)
        return result, already_seen

def get_all_connected_groups(graph):
    already_seen = set()
    result = []
    for node in graph:
        if node not in already_seen:
            connected_group, already_seen = get_connected_group(node, graph, already_seen)
            result.append(connected_group)
    return result
            
def connected_components(neighbors):
    seen = set()
    def component(node):
        nodes = set([node])
        while nodes:
            node = nodes.pop()
            seen.add(node)
            nodes |= neighbors[node] - seen
            yield node
    for node in neighbors:
        if node not in seen:
            yield component(node)

def n_components(graph):
    new_graph = {node: set(edge for edge in edges)
             for node, edges in graph.items()}
    components = []
    for component in connected_components(new_graph):
        c = set(component)
        components.append([edge for edges in graph.values()
                                for edge in edges
                                if c.intersection(edge)])

    return len(components)

bcs = [1,1,3,1,1,1,2,2,2,2,1]
fwd_umi = ['A','C','C','D','A','B','C','D','A','B','A']
rev_umi = ['A','B','D','D','B','B','C','D','A','B','A']

unique_tracker_either = {}


for line in range(len(bcs)):
    bc = bcs[line]
    umi1 = fwd_umi[line]
    umi2 = rev_umi[line]

    if bc not in unique_tracker_either.keys():
        unique_tracker_either[bc] = {}
        unique_tracker_either[bc]['F'+umi1] = set(['R'+umi2])
        unique_tracker_either[bc]['R'+umi2] = set(['F'+umi1])
    else:
        if 'F'+umi1 not in unique_tracker_either[bc].keys():
            unique_tracker_either[bc]['F'+umi1] = set(['R'+umi2])
        else:
            unique_tracker_either[bc]['F'+umi1].add('R'+umi2)
        if 'R'+umi2 not in unique_tracker_either[bc].keys():
            unique_tracker_either[bc]['R'+umi2] = set(['F'+umi1])
        else:
            unique_tracker_either[bc]['R'+umi2].add('F'+umi1)

print(unique_tracker_either)
for bc in unique_tracker_either.keys():
    components = n_components(unique_tracker_either[bc])

    print(components)
    # print(len(components))
    

        # f_out.write(bc + '\t' + str(len(components))+'\n')