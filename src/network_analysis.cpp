#include "network_analysis.h"

NetworkAnalysis::NetworkAnalysis(const std::string & filename_v,
                                char sep_v,
                                const string & filename_e,
                                char sep_e,
                                bool from_0)
{
    read_vertices(filename_v, sep_v);
    read_edges(filename_e, sep_e, from_0);
    print_network_info();
    std::cout << std::endl;
    init_network_complex();
    std::cout << std::endl;
    init_alpha_complex();
}


void NetworkAnalysis::read_vertices(const std::string & filename, char sep)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Could not open file: " << filename << std::endl;
        exit(-1);
    }
    string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        vector<double> point;
        double coord;
        while (ss >> coord) {
            point.push_back(coord);
            if (ss.peek() == sep)
                ss.ignore();
        }
        unsigned n_coord = point.size();
        Point p(n_coord, point.begin(), point.end());
        vertices.push_back(p);
    }
    file.close();
}


void NetworkAnalysis::read_edges(const std::string & filename, char sep, bool from_0)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Could not open file: " << filename << std::endl;
        exit(-1);
    }
    string line;
    unsigned u = from_0 ? 0 : 1;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        Edge edge; double index;
        ss >> index; edge[0]=index-u;
        if (ss.peek() == sep)
            ss.ignore();
        ss >> index; edge[1]=index-u;
        edges.push_back(edge);
    }
    file.close();
}


void NetworkAnalysis::init_network_complex()
{
    // Initialize the simplicial complex representing the network
    for(unsigned i = 0; i <edges.size(); ++i)
    {
        auto returnValue = network_simplex_tree.insert_simplex_and_subfaces(edges[i]);
        if(!returnValue.second)
            rej_edges_indexes.push_back(i);
    }
    // Display information about the network complex
    std::clog << "####################################" << std::endl;
    std::clog << "Network simplicial complex:" << std::endl;
    print_simplex_info(network_simplex_tree);
    std::clog << "* Number of rejected edges: " << rej_edges_indexes.size() << std::endl;
}


void NetworkAnalysis::init_alpha_complex()
{
    // Initialize an alpha complex from the list of points
    Gudhi::alpha_complex::Alpha_complex<Kernel> alpha_complex(vertices);
    if (alpha_complex.create_complex(alpha_simplex_tree))
    {
        // Display information about the alpha complex
        std::clog << "####################################" << std::endl;
        std::clog << "Alpha simplicial complex:" << std::endl;
        print_simplex_info(alpha_simplex_tree);
    }
    else
        std::cerr << "Error: alpha complex initialization fails" << std::endl;
}


void NetworkAnalysis::print_network_info() const
{
    std::cout << "####################################" << std::endl;
    std::cout << "Network formed by:\n* "<<
    vertices.size() << " vertices\n* " << edges.size() << " edges" << std::endl;
}

void NetworkAnalysis::print_betti_numbers() const
{
    std::clog << "The Betti numbers are:\n";
    unsigned i=0;
    for(unsigned b : betti_numbers)
        std::cout << "* b" << i++ << " = " << b << std::endl;
}

void NetworkAnalysis::print_max_voids_diameter(unsigned nmax) const
{
    unsigned n = std::min(nmax, static_cast<unsigned int>(voids_pers_int.size()));
    std::cout << "First " << n << " most persistent voids with extimated diameters:" << std::endl;
    std::cout << "(birth, death) -> diameter" << std::endl;
    for(unsigned i = 0; i<n; ++i)
        std::cout <<  "("<<voids_pers_int[i].first<<", "<<voids_pers_int[i].second<<") -> "<< max_voids_diameter[i] << std::endl;
}

void NetworkAnalysis::compute_network_homology()
{
    // By default, since the complex has dimension 1, only 0-dimensional homology would be computed.
    // Here we also want persistent homology to be computed for the maximal dimension in the complex (persistence_dim_max = true)
    Persistent_cohomology pcoh(network_simplex_tree, true);
    // Initialize the coefficient field Z/2Z for homology
    pcoh.init_coefficients(2);
    // Compute the persistence of the complex
    pcoh.compute_persistent_cohomology();
    // Compute Betti numbers
    betti_numbers = pcoh.betti_numbers();
    // Print them:
    print_betti_numbers();
}


void NetworkAnalysis::compute_alpha_persistent_homology(bool write_file)
{

        Persistent_cohomology pcoh(alpha_simplex_tree);
        // initializes the coefficient field for homology
        pcoh.init_coefficients(2);
        // Compute the persistence diagram of the complex
        pcoh.compute_persistent_cohomology();
        // calculate the maximum diameter of each void
        // get the persistence intervals for voids
        voids_pers_int = pcoh.intervals_in_dimension(alpha_simplex_tree.dimension()-1);
        // sort the persistence intervals by decreasing life time For
        std::sort(begin(voids_pers_int),
                  end(voids_pers_int),
                  []( auto const& p1, auto const& p2){return (p1.second-p1.first) > (p2.second-p2.first);} );
        // resize the diameters vector
        max_voids_diameter.resize(voids_pers_int.size());
        // store diameters
        std::transform(cbegin(voids_pers_int),
                       cend(voids_pers_int),
                       begin(max_voids_diameter),
                       []( auto const& p){return 2 * std::sqrt(p.second);} );

        // write persistence in a file
        if(write_file)
        {
            std::stringstream ss;
            ss << "../../output/alpha_persistence";
            std::ofstream out(ss.str().c_str());
            pcoh.output_diagram(out);
            out.close();
        }
}


void print_simplex_info(Simplex_tree & st)
{
    std::clog << "* Dimension " << st.dimension() << std::endl;
    std::clog << "* Contains " << st.num_simplices() << " simplices" << std::endl;
    std::clog << "* Contains " << st.num_vertices() << " vertices" << std::endl;
}


void print_range_simplices(Simplex_tree & st, unsigned start_index, unsigned length)
{
    unsigned n = st.num_simplices();
    if(start_index < 0 or start_index + length > n-1)
    {
        std::cerr << "Error: invalid range, index out of bound" << std::endl;
        return;
    }
    std::clog << "Iterator on simplices from "<<start_index<<" to " << start_index+length-1<<", with [filtration value]:"<<std::endl;
    auto it = st.filtration_simplex_range().cbegin() + start_index;
    unsigned i=0;
    for(; i < length; ++it,++i)
    {
        std::clog << "   " << "[" << st.filtration(*it) << "] ";
        for (const auto & vertex : st.simplex_vertex_range(*it))
            std::clog << "(" << vertex << ")";
        std::clog << std::endl;
    }
}
