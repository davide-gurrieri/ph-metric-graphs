import sys
sys.path.append('../src')
from network import *
from utils import *
import matplotlib.pyplot as plt

net1 = Network()
net1.read_vertices(name = "../../data/complex_network/VertexCoordinates.txt",
                 sep = ' ')
net1.read_edges(name = "../../data/complex_network/EdgeConnectivity.txt",
                 sep = ' ')

alpha_complex = gd.AlphaComplex(points = net1.vertices, precision="safe")
alpha_simplex_tree = alpha_complex.create_simplex_tree()

dim = alpha_simplex_tree.dimension()

alpha_barcode = alpha_simplex_tree.persistence()

for i in range(dim): # 0 1 ... dim-1
    barcode_i = alpha_simplex_tree.persistence_intervals_in_dimension(i)
    col=plt.cm.Set1.colors[i:]
    gd.plot_persistence_barcode(barcode_i, max_intervals=1000, colormap=col, alpha=1)
    plt.savefig(f'../../output/plot/barcode_dim_{i}.svg')
    plt.savefig(f'../../output/plot/barcode_dim_{i}.png')

gd.plot_persistence_diagram(alpha_barcode,alpha=1,legend=True)
plt.savefig(f'../../output/plot/persistence_diagram.png')
plt.savefig(f'../../output/plot/persistence_diagram.svg')
