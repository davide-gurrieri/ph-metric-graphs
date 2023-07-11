import gudhi as gd
from network import *


class NetworkHomology:
    def __init__(self, network):
        self.network = network

        self.simplex_tree = gd.SimplexTree()
        # save 1-simplices' indexes that are not stored in the graph_simplex_tree
        # (many copies or degenerate edge)
        self.rejected_indexes = []
        rejected_edges = []
        i = 0
        for edge in network.edges:
            inserted = self.simplex_tree.insert(edge)
            if not inserted:
                self.rejected_indexes.append(i)
                rejected_edges.append(edge)
            i+=1

        self.n_rejected = len(self.rejected_indexes)
        self.rejected_edges = np.array(rejected_edges)
        self.dim = self.simplex_tree.dimension()
        self.n_simplices = self.simplex_tree.num_simplices()
        # compute homology
        self.simplex_tree.set_dimension(2)
        self.simplex_tree.compute_persistence()
        self.betti_numbers = self.simplex_tree.betti_numbers()

    def get_simplices_dim(self, dim):
        return np.array([s[0] for s in self.simplex_tree.get_skeleton(dim) if len(s[0]) == dim+1])

    def __str__(self):
        n_0 = self.simplex_tree.num_vertices()
        n_1 = self.get_simplices_dim(1).shape[0]
        return (f"Simplicial complex of dimension {self.dim}:\n"
                f"Number of simplices: {self.n_simplices}\n"
                f"      {n_0} of dimension 0\n"
                f"      {n_1} of dimension 1\n"
                f"{self.n_rejected} edges were discarded\n"
                f"0th Betti number: {self.betti_numbers[0]}\n"
                f"1st Betti number: {self.betti_numbers[1]}\n"
                f"Interpretation: the network is composed of {self.betti_numbers[0]} connected components and has {self.betti_numbers[1]} loops")
