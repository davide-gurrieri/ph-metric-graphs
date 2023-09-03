#include "network_analysis.h"

void Filtration::compute_persistent_cohomology(bool persistence_dim_max) {
  persistent_cohomology = std::make_unique<Persistent_cohomology>(
      simplex_tree, persistence_dim_max);
  // Initialize the coefficient field Z/2Z for homology
  persistent_cohomology->init_coefficients(2);
  // Compute the persistence of the complex
  chrono.start();
  persistent_cohomology->compute_persistent_cohomology();
  chrono.stop();
}

void Filtration::save_persistence(const string &filename) const {
  std::ofstream out(filename.c_str());
  persistent_cohomology->output_diagram(out);
  out.close();
}

void Filtration::print_complex_info() {
  std::cout << "Simplicial complex information:" << std::endl;
  std::cout << "* Dimension " << simplex_tree.dimension() << std::endl;
  std::cout << "* Contains " << simplex_tree.num_simplices() << " simplices"
            << std::endl;
  std::cout << "* Contains " << simplex_tree.num_vertices() << " vertices"
            << std::endl;
}

void Filtration::print_range_simplices(unsigned start_index, unsigned length) {
  unsigned n = simplex_tree.num_simplices();
  if (start_index < 0 or start_index + length > n - 1) {
    std::cerr << "Error: invalid range, index out of bound" << std::endl;
    return;
  }
  std::clog << "Simplices from " << start_index << " to "
            << start_index + length - 1
            << ", with [filtration value]:" << std::endl;
  auto it = simplex_tree.filtration_simplex_range().cbegin() + start_index;
  unsigned i = 0;
  for (; i < length; ++it, ++i) {
    std::clog << "   "
              << "[" << simplex_tree.filtration(*it) << "] ";
    for (const auto &vertex : simplex_tree.simplex_vertex_range(*it))
      std::clog << "(" << vertex << ")";
    std::clog << std::endl;
  }
}

void Filtration::print_betti_numbers() const {
  std::cout << "Betti numbers of the complete simplicial complex are:\n";
  unsigned i = 0;
  for (unsigned b : persistent_cohomology->betti_numbers())
    std::cout << "* b" << i++ << " = " << b << std::endl;
}

void Filtration::print_elapsed_time() const {
  std::cout << "Elapsed time for persistent homology computation: "
            << chrono.wallTime() * 0.001 << " milliseconds" << std::endl;
}

// ##########################################################################

AlphaFiltration::AlphaFiltration(const vector<Point<CoordType>> &points) {
  // Initialize an alpha complex from the list of points
  vector<AlphaPoint> alpha_points(points.size());
  std::transform(cbegin(points), cend(points), begin(alpha_points),
                 [](auto const &p) {
                   return AlphaPoint(p.dimension(), p.coordinates.begin(),
                                     p.coordinates.end());
                 });
  Alpha_complex alpha_complex(alpha_points);
  if (!alpha_complex.create_complex(simplex_tree)) {
    std::cerr << "Error: alpha complex initialization fails" << std::endl;
    exit(1);
  }
}

void AlphaFiltration::compute_voids_diameter() {
  // calculate the maximum diameter of each void
  // get the persistence intervals for voids
  voids_persistence = persistent_cohomology->intervals_in_dimension(
      simplex_tree.dimension() - 1);
  // sort the persistence intervals by decreasing life time
  std::sort(begin(voids_persistence), end(voids_persistence),
            [](auto const &p1, auto const &p2) {
              return (p1.second - p1.first) > (p2.second - p2.first);
            });
  // resize the diameters vector
  voids_diameter.resize(voids_persistence.size());
  // store diameters
  std::transform(cbegin(voids_persistence), cend(voids_persistence),
                 begin(voids_diameter),
                 [](auto const &p) { return 2 * std::sqrt(p.second); });
}

void AlphaFiltration::print_diameters(unsigned nmax) const {
  unsigned n = std::min(nmax, static_cast<unsigned int>(voids_diameter.size()));
  std::cout << "First " << n
            << " most persistent voids with extimated diameters:" << std::endl;
  std::cout << "(birth, death); sqrt(persistence); diameter" << std::endl;
  for (unsigned i = 0; i < n; ++i)
    std::cout << "(" << voids_persistence[i].first << ", "
              << voids_persistence[i].second << "); "
              << std::sqrt(voids_persistence[i].second -
                           voids_persistence[i].first)
              << "; " << voids_diameter[i] << std::endl;
}

void AlphaFiltration::make_analysis() {
  compute_persistent_cohomology(false);
  compute_voids_diameter();
}

void AlphaFiltration::print_analysis() {
  std::cout << "#################################################" << std::endl;
  std::cout << "ALPHA FILTRATION ANALYSIS" << std::endl;
  std::cout << "#################################################" << std::endl;
  std::cout << std::endl;
  print_complex_info();
  std::cout << std::endl;
  print_betti_numbers();
  std::cout << std::endl;
  print_diameters(10);
  std::cout << std::endl;
  print_elapsed_time();
  std::cout << std::endl;
}

// ##########################################################################

RadialFiltration::RadialFiltration(const vector<Point<CoordType>> &points,
                                   const vector<Edge> &edges) {
  if (points.empty()) {
    std::cerr << "Error: empty point cloud" << std::endl;
    exit(1);
  }
  Point<CoordType> center(points[0].dimension());
  for (const auto &p : points)
    center = center + p;
  center = center / points.size();
  double filtration_value;
  for (int i = 0; i < points.size(); ++i) {
    filtration_value = (points[i] - center).norm();
    if (filtration_value > max_radius)
      max_radius = filtration_value;

    simplex_tree.insert_simplex({i}, filtration_value);
  }
  for (const auto &edge : edges) {
    const double filtration_value = std::max((points[edge[0]] - center).norm(),
                                             (points[edge[1]] - center).norm());
    simplex_tree.insert_simplex_and_subfaces(edge, filtration_value);
  }
  n_edges = simplex_tree.num_simplices() - simplex_tree.num_vertices();
}

void RadialFiltration::compute_descriptors() {
  // compute the descriptors
  std::vector<std::pair<double, double>> persistence_0 =
      persistent_cohomology->intervals_in_dimension(0);

  tortuosity_descriptor =
      std::count_if(persistence_0.begin(), persistence_0.end(), [this](auto p) {
        return (p.second - p.first) <= 0.1 * max_radius;
      });
  tortuosity_descriptor /= n_edges;
  loops_descriptor =
      persistent_cohomology->betti_number(1) / static_cast<double>(n_edges);
}

void RadialFiltration::make_analysis() {
  compute_persistent_cohomology(true);
  compute_descriptors();
}

void RadialFiltration::print_analysis() {
  std::cout << "#################################################" << std::endl;
  std::cout << "RADIAL FILTRATION ANALYSIS" << std::endl;
  std::cout << "#################################################" << std::endl;
  std::cout << std::endl;
  print_complex_info();
  std::cout << std::endl;
  print_betti_numbers();
  std::cout << std::endl;
  std::cout << "Tortuosity descriptor: " << tortuosity_descriptor << std::endl;
  std::cout << "Loops descriptor: " << loops_descriptor << std::endl;
  std::cout << std::endl;
  print_elapsed_time();
  std::cout << std::endl;
}

// ##########################################################################

void NetworkData::print_network_info() const {
  std::cout << "#################################################" << std::endl;
  std::cout << std::endl;
  std::cout << "Input vascular network composed of:\n* " << vertices.size()
            << " vertices\n* " << edges.size() << " edges" << std::endl;
  std::cout << std::endl;
}

void NetworkData::parse(const std::string &filename_v, char sep_v,
                        const string &filename_e, char sep_e, bool from_0) {
  read_vertices(filename_v, sep_v);
  read_edges(filename_e, sep_e, from_0);
}

void NetworkData::read_vertices(const std::string &filename, char sep) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Could not open file: " << filename << std::endl;
    exit(-1);
  }
  string line;
  while (std::getline(file, line)) {
    std::stringstream ss(line);
    vector<CoordType> point;
    CoordType coord;
    while (ss >> coord) {
      point.push_back(coord);
      if (ss.peek() == sep)
        ss.ignore();
    }
    vertices.push_back(Point<CoordType>(point));
  }
  file.close();
}

void NetworkData::read_edges(const std::string &filename, char sep,
                             bool from_0) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Could not open file: " << filename << std::endl;
    exit(-1);
  }
  string line;
  unsigned u = from_0 ? 0 : 1;
  while (std::getline(file, line)) {
    std::stringstream ss(line);
    Edge edge;
    double index;
    ss >> index;
    edge[0] = index - u;
    if (ss.peek() == sep)
      ss.ignore();
    ss >> index;
    edge[1] = index - u;
    edges.push_back(edge);
  }
  file.close();
}

// ##########################################################################

NetworkAnalysis::NetworkAnalysis(const NetworkData *data_) : data(data_) {
  filtrations.push_back(
      std::make_unique<RadialFiltration>(data->vertices, data->edges));

  filtrations.push_back(std::make_unique<AlphaFiltration>(data->vertices));
}

void NetworkAnalysis::analyse() const {
  for (auto &filtration : filtrations)
    filtration->make_analysis();
}

void NetworkAnalysis::print_global_analysis() const {
  for (auto &filtration : filtrations)
    filtration->print_analysis();
}

void NetworkAnalysis::save_persistence(const std::string &name) const {
  std::vector<std::string> names = {"radial", "alpha"};
  unsigned i = 0;
  for (const auto &filtration : filtrations)
    filtration->save_persistence("../../output/" + name + "_persistence_" +
                                 names[i++]);
}