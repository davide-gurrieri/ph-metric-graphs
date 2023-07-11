import sys
sys.path.append('../src')
from network import *
from utils import *

net1 = Network()
net1.read_vertices(name = "../../data/complex_network/VertexCoordinates.txt",
                 sep = ' ')
net1.read_edges(name = "../../data/complex_network/EdgeConnectivity.txt",
                 sep = ' ')

alpha_complex = gd.AlphaComplex(points = net1.vertices, precision="safe")
alpha_simplex_tree = alpha_complex.create_simplex_tree()

_ = plot_simplicial_complex(
    net1.vertices,
    alpha_simplex_tree,
    filtration_value=20000,
    points_thickness = 0.7,
    points_colour = 'midnightblue',
    segments_thickness = 0.7,
    segments_colour = 'midnightblue',
    other_color = 'salmon',
    show = True
)
