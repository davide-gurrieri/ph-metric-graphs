/**
 * @file network_analysis.h
 * @author Davide Gurrieri (davide.gurrieri@mail.polimi.it)
 * @brief This file contains the declaration of the NetworkAnalysis class.
 * @version 0.1
 * @date 2023-08-29
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <gudhi/Alpha_complex.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Simplex_tree.h>

#include <CGAL/Epeck_d.h>

using std::string;
using std::vector;

using Kernel = CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>;

using Point = Kernel::Point_d;
using Edge = std::array<unsigned, 2>;

using Simplex_tree = Gudhi::Simplex_tree<>;
using Filtration_value = Simplex_tree::Filtration_value;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology =
    Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;
using Alpha_complex = Gudhi::alpha_complex::Alpha_complex<Kernel>;

class NetworkAnalysis {
private:
  vector<Point> vertices;
  vector<Edge> edges;
  vector<double> edges_length;
  vector<unsigned> rej_edges_indexes;

public:
  Simplex_tree network_simplex_tree;
  Simplex_tree alpha_simplex_tree;
  vector<int> betti_numbers;
  vector<double> max_voids_diameter;
  vector<std::pair<Filtration_value, Filtration_value>> voids_pers_int;

  /**
   * @brief Default constructor
   */
  NetworkAnalysis() = default;

  /**
   * @brief User defined constructor that reads data from files
   *
   * The constructor calls two member functions of the class: read_vertices
   * and read_edges passing them proper parameters. Then, calls
   * init_network_complex and init_alpha_complex to initialize the network and
   * alpha simplex trees, respectively.
   *
   * @param filename_v the name of a file that contains the vertex data for the
   * network (including path)
   * @param sep_v the separator character used in the vertex data file (like ','
   * or ' ')
   * @param filename_e the name of a file that contains the edge data for the
   * network. (including path)
   * @param sep_e the separator character used in the edge data file (like ','
   * or ' ')
   * @param from_0 a boolean flag that indicates whether the vertex indices in
   * the edge data file are 0-based (i.e., start from 0) or 1-based (i.e., start
   * from 1)
   */
  NetworkAnalysis(const std::string &filename_v, char sep_v,
                  const string &filename_e, char sep_e, bool from_0 = true);

  /**
   * @brief Computes the homology of the network (number of connected components
   * and number of loops)
   *
   * Then initializes the coefficient field for homology to be
   * $\mathbb{Z}/2\mathbb{Z}$, and computes the persistent cohomology of the
   * network simplex tree. The Betti numbers of the network are printed to the
   * console.
   *
   */
  void compute_network_homology();

  /**
   * @brief Computes the alpha persistent homology for a given alpha simplex
   * tree and stores voids ph informations and an exstimate of their maximum
   * diameter.
   *
   * The function uses the alpha simplex tree to compute the alpha persistent
   * homology. PH is then used to calculate an extimate of the maximum diameter
   * of each void in the complex. The filtration parameters corresponding to
   * birth and death of each voids and their diameters are stored in the class
   * members voids_pers_int and max_voids_diameter, respectively. If the
   * write_file parameter is set to true, the persistence informations are
   * written to a file named "alpha_persistence" in the "../output" directory.
   * This file could be usefull to plot persistence barcodes and diagrams.
   *
   * @param write_file A boolean value indicating whether to write the
   * persistence diagram to a file.
   * @return void
   */
  void compute_alpha_persistent_homology(bool write_file);

  /**
   * @brief getter for the vertices vector of the network
   * @param void
   * @return vector<CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>::Point_d>
   */
  vector<Point> get_vertices() const { return vertices; }

  /**
   * @brief getter for the edges vector of the network
   * @param void
   * @return vector<array<unsigned,2>>
   */
  vector<Edge> get_edges() const { return edges; }

  /**
   * @brief getter for the vector of the edges' indexes not inserted in the
   * network simplex tree
   *
   * Edges that degenerate into a vertices or alredy in the simplex are not
   * inserted.
   *
   * @param void
   * @return vector<unsigned>
   */
  vector<unsigned> get_rej_edges_indexes() const { return rej_edges_indexes; }

  void print_network_info() const;
  void print_betti_numbers() const;
  void print_max_voids_diameter(unsigned nmax) const;

private:
  /**
   * @brief Reads vertex data from a file
   *
   * Reads vertex data from a file using std::getline,
   * extracts the coordinate data using a std::stringstream, and stores it in a
   * Point object (i.e. CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>::Point_d). If
   * the next character in the stream is equal to `sep`, it is ignored using the
   * ignore() function. Points are stored in a class member vector `vertices`.
   * @param filename the name of a file that contains the vertex data for the
   * network (including path)
   * @param sep the separator character used in the vertex data file (like ','
   * or ' ')
   */
  void read_vertices(const string &filename, char sep);

  /**
   * @brief Reads edge data from a file
   *
   * Similar to read_vertices. In this case the vertices' indices of each edge
   * are stored in an Edge (i.e std::array<unsigned, 2>) and are decremented
   * by 1 if 1-based.
   * Edges are stored in a class member vector `edges`.
   *
   * @param filename the name of a file that contains the edge data for the
   * network. (including path)
   * @param sep the separator character used in the edge data file (like ',' or
   * ' ')
   * @param from_0 a boolean flag that indicates whether the vertex indices in
   * the edge data file are 0-based (i.e., start from 0) or 1-based (i.e., start
   * from 1)
   */
  void read_edges(const string &filename, char sep, bool from_0);

  /**
   * @brief initialize the simplicial complex representing the network
   *
   * Iterates over each edge in the edges vector, inserts it into the network
   * simplex tree data structure, and keeps track of any rejected edges by
   * pushing their indexes onto the `rej_edges_indexes` vector.
   *
   * @param void
   * @return void
   */
  void init_network_complex();

  /**
   * @brief Initializes an alpha complex from a list of points and displays
   * information about the complex.
   *
   * Initializes an alpha complex from the vector of Points `vertices`, using
   * the Gudhi library. If the alpha complex is successfully created, displays
   * informations about the complex. Otherwise prints an error message.
   *
   * @param void
   * @return void
   */
  void init_alpha_complex();
};

/**
 * @brief Prints a range of simplices from a simplex tree, with their filtration
 * values
 *
 * Prints to the console the simplices in the range [start_index,
 * start_index+length-1] of the simplex tree st, along with their filtration
 * values. The function returns an error message and does not print anything if
 * the range is invalid.
 *
 * @param st The simplex tree to print from
 * @param start_index The index of the first simplex to print
 * @param length The number of simplices to print
 * @return void
 */
void print_range_simplices(Simplex_tree &st, unsigned start_index,
                           unsigned length);

/**
 * @brief Prints useful informations about the simplex tree passed as input
 * @param st The simplex tree to print from
 * @return void
 */
void print_simplex_info(Simplex_tree &st);

/*
Using the default CGAL::Epeck_d makes the construction safe.
If you pass exact=true to create_complex, the filtration values are the exact
ones converted to the filtration value type of the simplicial complex. This can
be very slow. If you pass exact=false (the default), the filtration values are
only guaranteed to have a small multiplicative error compared to the exact
value.

Using CGAL::Epick_d makes the computations slightly faster, and the
combinatorics are still exact, but the computation of filtration values can
exceptionally be arbitrarily bad. In all cases, it is still guaranteed that the
output is a valid filtration (faces have a filtration value no larger than their
cofaces).
*/
