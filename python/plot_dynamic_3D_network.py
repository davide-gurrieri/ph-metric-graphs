from network import *

net1 = Network()
net1.read_vertices(name = "../data/complex_network/VertexCoordinates.txt", sep = ' ')
net1.read_edges(name = "../data/complex_network/EdgeConnectivity.txt", sep = ' ')

_ = net1.plot(segments_thickness = 1,
            segments_color = 'black',
            vertices_thickness = 0.5,
            vertices_color = 'red',
            show = True)
