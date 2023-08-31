from network import *

# Synthetic networks
net = Network()
for i in range(0,7):
    net.read_vertices(name = f"../data/synthetic_network/VerticesCoordinates_negSeed0.{i}.txt", sep = ',')
    net.read_edges(name = f"../data/synthetic_network/Connectivity_negSeed0.{i}.txt", sep = ',')
    _ = net.plot_static(segments_thickness = 1.5,
            segments_color = 'black',
            vertices_thickness = 4,
            vertices_color = 'red',
            show = False)
            
    plt.savefig(f"../output/plot/syntetic_network_{i}.pdf")
    plt.close()
    
# Real network
net1 = Network()
net1.read_vertices(name = "../data/complex_network/VertexCoordinates.txt", sep = ' ')
net1.read_edges(name = "../data/complex_network/EdgeConnectivity.txt", sep = ' ')

_ = net1.plot_static(segments_thickness = 0.1,
            segments_color = 'black',
            vertices_thickness = 0.005,
            vertices_color = 'red',
            show = False)

plt.savefig(f"../output/plot/real_network.pdf")
plt.close()