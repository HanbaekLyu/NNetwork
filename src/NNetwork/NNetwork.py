import numpy as np

"""
Neighborhood Network package
Simple and weighted network objects with built-in neighborhood list for scalable random walk and MCMC sampling
Author: Josh Vendrow and Hanbaek Lyu
"""

class NNetwork():

    def __init__(self):
        self.edges = []
        # self.sorted_left = None
        # self.sorted_right = None
        self.neighb = {}
        self.vertices = []
        self.vertices_set = set()
        self.number_nodes = 0

    """ 
    def set_edges(self, edges):

        self.sort_left(edges)
        self.sort_right(edges)
        self.edges = edges
        self.set_neighbors()
    """

    def add_edges(self, edges):
        """Given an edgelist, add each edge in the edgelist to the network"""
        i = 0
        for edge in edges:
            self.add_edge(edge, update=False)
            # print('Loading %i th edge out of %i edges' % (i, len(edges)))
            # i += 1

        # self.node = list(self.neighb.keys())

    def add_edge(self, edge, update=True):
        """Given an edge, add this edge to the Network"""

        u, v = edge
        self.edges.append(edge)
        if u not in self.neighb:
            self.neighb[u] = set({v})
            self.vertices.append(u)
            self.vertices_set.add(u)
            self.number_nodes += 1

        else:
            self.neighb[u].add(v)

        if v not in self.neighb:
            self.neighb[v] = set({u})
            self.vertices.append(v)
            self.vertices_set.add(v)
            self.number_nodes += 1

        else:
            self.neighb[v].add(u)

            
    def delete_edge(self, edge,):
        """Delete edge from edgelist"""

        u, v = edge

        if u in self.neighb and v in self.neighb[u]:
            self.neighb[u].remove(v)

        if v in self.neighb and u in self.neighb[v]:
            self.neighb[v].remove(u)

        for i in range(len(self.edges)):
            if self.edges[i] == edge:
                self.edges.pop(i)
                
                
    def add_nodes(self, nodes):
        """Given a list of nodes, adds the nodes to the Network"""

        for node in nodes:
            self.add_node(node)

    def add_node(self, node):
        """Given a single node, adds the node to the Network"""

        if node not in self.vertices:
            self.vertices.append(node)
            self.neighb[node] = set()
            self.number_nodes += 1
            self.vertices_set.add(node)

    def set_neighbors(self):

        for edge in self.edges:
            self.add_edge(edge)

    def get_edges(self):
        set_edgelist = []
        for x in self.vertices:
            if x in self.neighb:
                for y in self.neighb[x]:
                    set_edgelist.append([x, y])
        self.edges = set_edgelist
        return self.edges

        # self.node = list(self.neighb.keys())

    def intersection(self, Network_other):
        """
        Given another network, returns all edges found in both networks.
        """
        common_nodeset = self.vertices_set.intersection(Network_other.vertices_set)
        common_edgelist = []

        for x in common_nodeset:
            for y in self.neighb[x].intersection(Network_other.neighbors(x)):
                common_edgelist.append([x, y])

        return common_edgelist

    def neighbors(self, node):
        """
        Given a node, returns all the neighbors of the node.
        """
        return self.neighb[node]

    def has_edge(self, u, v):
        """
        Given two nodes, returns true of these is an edge between then, false otherwise.
        """
        try:
            return v in self.neighb[u]
        except KeyError:
            return False

    def nodes(self, is_set=False):
        return self.vertices if not is_set else self.vertices_set

    def num_nodes(self):
        return self.number_nodes

    def edges(self):
        return self.edges


class Wtd_NNetwork():
    '''
    Class for weighted and undirectional network
    Can store additional tensor information in colored_edge_weights -- list of nonnegative numbers
    Underlying NNetwork serves as the skeleton for RWs and motif sampling
    colored edges = additional information on the (already weighted) edges
    Think of edges = pixel location and colored edges = color of that pixel
    '''

    def __init__(self):
        self.edges = []
        self.edge_weights = {}
        self.neighb = {}
        self.vertices = []
        self.vertices_set = set()
        self.number_nodes = 0
        self.colored_edge_weights = {}
        self.color_dim = 0

    def add_edges(self, edges, edge_weight=1, increment_weights=True):
        """Given an edgelist, add each edge in the edgelist to the network"""
        for edge in edges:
            self.add_edge(edge, weight=edge_weight, increment_weights=increment_weights)
            # self.edge_weights.update({str(edge): default_edge_weight})

    def add_edge(self, edge, weight=1, increment_weights=False):
        """Given an edge, add this edge to the Network"""
        u, v = edge
        # self.edges.append(edge)
        if not increment_weights or not self.has_edge(u,v):
            wt = weight
        else:
            wt = self.get_edge_weight(u, v) + weight

        self.edge_weights.update({str([str(u), str(v)]): wt})
        self.edge_weights.update({str([str(v), str(u)]): wt})

        if u not in self.neighb:
            self.neighb[u] = {v}
            self.vertices.append(u)
            self.vertices_set.add(u)
            self.number_nodes += 1
        else:
            self.neighb[u].add(v)

        if v not in self.neighb:
            self.neighb[v] = {u}
            self.vertices.append(v)
            self.vertices_set.add(v)
            self.number_nodes += 1
        else:
            self.neighb[v].add(u)

    def add_nodes(self, nodes):
        """Given a list of nodes, adds the nodes to the Network"""
        for node in nodes:
            self.add_node(node)

    def add_node(self, node):
        """Given a single node, adds the node to the Network"""
        if node not in self.vertices:
            self.vertices.append(node)
            self.neighb[node] = set()
            self.number_nodes += 1
            self.vertices_set.add(node)

    def set_neighbors(self):
        for edge in self.edges:
            self.add_edge(edge)

    def intersection(self, Network_other):
        """
        Given another network, returns all edges found in both networks.
        """
        edges_other = set(Network_other.edges())
        edges = set(self.get_edges())
        intersection = edges.intersection(edges_other)
        return list(intersection)

    def neighbors(self, node):
        """
        Given a node, returns all the neighbors of the node.
        """
        return self.neighb[node]

    def has_edge(self, u, v):
        """
        Given two nodes, returns true of these is an edge between then, false otherwise.
        """
        try:
            return v in self.neighb[u]
        except KeyError:
            return False

    def nodes(self, is_set=False):
        return self.vertices if not is_set else self.vertices_set

    def num_nodes(self):
        return self.number_nodes

    def get_edge_weight(self, u, v):
        return self.edge_weights.get(str([str(u), str(v)]))

    def get_edges(self):
        set_edgelist = []
        for x in self.vertices:
            if x in self.neighb:  ### x has a neighbor
                for y in self.neighb[x]:
                    set_edgelist.append([x, y])
        self.edges = set_edgelist
        return self.edges

    def get_wtd_edgelist(self):
        edgelist = self.get_edges()
        wtd_edgelist = []
        for edge in edgelist:
            u, v = edge
            wtd_edgelist.append([u, v, abs(self.get_edge_weight(u, v))])

        return wtd_edgelist

    def get_abs_wtd_edgelist(self):
        edgelist = self.get_edges()
        wtd_edgelist = []
        for edge in edgelist:
            u, v = edge
            wtd_edgelist.append([u, v, self.get_edge_weight(u, v)])

        return wtd_edgelist

    def save_wtd_edgelist(self, default_folder = 'Temp_save_graphs', default_name='temp_wtd_edgelist'):
        edgelist = self.get_edges()
        wtd_edgelist = []
        for edge in edgelist:
            u, v = edge
            wtd_edgelist.append([u, v, self.get_edge_weight(u, v)])

        ### Todo: Use pickle later
        with open(default_folder + '/' + default_name + '.txt', "w") as file:
            file.write(str(wtd_edgelist))
        return wtd_edgelist

    def add_wtd_edges(self, edges, increment_weights=False):
        """Given an edgelist, add each edge in the edgelist to the network"""
        for wtd_edge in edges:
            self.add_edge(wtd_edge[0:2], weight=wtd_edge[-1], increment_weights=increment_weights)

    def load_add_wtd_edges(self, path, increment_weights=False):
        with open(path, "r") as file:
            wtd_edgelist = eval(file.readline())
        self.add_wtd_edges(edges=wtd_edgelist, increment_weights=increment_weights)

    def clear_edges(self):
        node_set = self.vertices
        self.__init__()
        self.add_nodes(node_set)

    def threshold2simple(self, threshold=0.5):
        ### Threshold edges of the given weighted network and return a simple graph
        G_simple = NNetwork()
        for edge in self.get_edges():
            u, v = edge
            edge_weight = self.get_edge_weight(u, v)
            if edge_weight > threshold:
                G_simple.add_edge([u, v])
        return G_simple

    def add_colored_edges(self, colored_edges):
        """Given a list of colored edges, add each colored edge to the network"""
        for edge in colored_edges:
            self.add_colored_edge(edge)

    def add_colored_edge(self, colored_edge):
        """Given a colored edge, [u, v, a1, a2, ...], add this to the Network"""
        u, v = colored_edge[0:2]

        self.colored_edge_weights.update({str([str(u), str(v)]): colored_edge[2:]})
        self.colored_edge_weights.update({str([str(v), str(u)]): colored_edge[2:]})

    def set_clrd_edges_signs(self):
        """
        colored edge weight of and edge with weight w = [+(w), -(w)] (split positive and negative parts as a vector)
        """
        edgelist = self.get_edges()

        for edge in edgelist:
            u, v = edge
            w = self.get_edge_weight(u, v)
            colored_edge = [u, v, w*(w>=0), -w*(w<0)]
            self.add_colored_edge(colored_edge)

        u, v = self.get_edges()[0]
        self.color_dim = len(self.get_colored_edge_weight(u, v))

    def get_colored_edge_weight(self, u, v):
        return self.colored_edge_weights.get(str([str(u), str(v)]))

