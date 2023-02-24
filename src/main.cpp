#include "network_analysis.h"
#include "GetPot"
#include "chrono.hpp"

int main(int argc, char* argv[])
{
    GetPot command_line(argc, argv);
    const string file_name = command_line("file", "../options1.txt");
    GetPot datafile(file_name.c_str());

    string vertices_file_name = datafile("vertices_file_name", "../");
    string edges_file_name = datafile("edges_file_name", "../");
    char vertex_sep = datafile("vertex_sep", ',');
    char edge_sep = datafile("edge_sep", ',');
    bool indexes_from_0 = datafile("indexes_from_0", 0);


    Timings::Chrono chrono;

    // construct a NetworkAnalysis object from files
    NetworkAnalysis net = NetworkAnalysis(vertices_file_name,
                                          vertex_sep,
                                          edges_file_name,
                                          edge_sep,
                                          indexes_from_0);


    std::cout<<std::endl;
    std::cout << "####################################" << std::endl;
    std::cout << "First 5 simplices of the network simplicial complex:" << std::endl;
    print_range_simplices(net.network_simplex_tree, 0, 5);
    std::cout<<std::endl;
    std::cout << "####################################" << std::endl;

    chrono.start();
    net.compute_network_homology();
    chrono.stop();
    std::cout << "Elapsed time for network homology computation: " << chrono.wallTime() << " microsec"  << std::endl;

    chrono.start();
    net.compute_alpha_persistent_homology(true);
    chrono.stop();
    std::cout << "Elapsed time for alpha persistent homology computation: " << chrono.wallTime() << " microsec" << std::endl;

    std::cout<<std::endl;
    std::cout << "####################################" << std::endl;
    net.print_max_voids_diameter(10);

    return 0;
}
