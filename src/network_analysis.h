/**
 * @file network_analysis.h
 * @author Davide Gurrieri (davide.gurrieri@mail.polimi.it)
 * @brief This file contains classes declarations needed for the
 * topology analysis of vascular networks.
 * @date 2023-08-29
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include "chrono.hpp"
#include "geometry.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
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
using CoordType = double;

// ##########################################################################

/**
 * @brief Abstract class that represent a generic filtration
 *
 */
class Filtration {
public:
  /// Represents a filtered simplicial complex
  using SimplexTree = Gudhi::Simplex_tree<>;
  using Filtration_value = SimplexTree::Filtration_value;
  /// Homology coefficents
  using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
  using Persistent_cohomology =
      Gudhi::persistent_cohomology::Persistent_cohomology<SimplexTree,
                                                          Field_Zp>;

protected:
  SimplexTree simplex_tree;
  std::unique_ptr<Persistent_cohomology> persistent_cohomology = nullptr;
  Chrono chrono;

public:
  /**
   * @brief Computes the persistent homology of the filtered simplicial complex
   *
   * @param persistence_dim_max Boolean flag that indicates whether to compute
   * the persistent homology also for the homology class of the dimension of the
   * complex ($N$). If false, only the persistent homology of classes of
   * dimension
   * $\leq n$ are computed.
   */
  void compute_persistent_cohomology(bool persistence_dim_max);

  /**
   * @brief Saves the persistence information to a file
   *
   * @param filename Name of the file (including path) where to save the
   * persistence diagram
   */
  void save_persistence(const string &filename) const;

  /**
   * @brief Prints information about the simplicial complex
   *
   * It cannot be const because SimplexTree::dimension() and
   * SimplexTree::num_simplices() are (incorrectly) not const
   */
  void print_complex_info();

  /**
   * @brief Debug function to print simplices and filtration values in a range
   *
   * @param start_index Index of the first simplex to print
   * @param length Number of simplices to print
   */
  void print_range_simplices(unsigned start_index, unsigned length);

  /**
   * @brief Prints the Betti numbers of the simplicial complex
   *
   * It cannot be const due to SimplexTree methods
   */
  void print_betti_numbers() const;

  /**
   * @brief prints the time taken to calculate the persistent homology
   *
   */
  void print_elapsed_time() const;

  /**
   * @brief Pure virtual function that performs the topological analysis
   *
   */
  virtual void make_analysis() = 0;

  /**
   * @brief Pure virtual function that prints the results of the topological
   * analysis
   *
   */
  virtual void print_analysis() = 0;

  virtual ~Filtration() = default;
};

// ##########################################################################

/**
 * @brief Derived class from Filtration. It represents an alpha filtration
 *
 */
class AlphaFiltration : public Filtration {
public:
  using Kernel = CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>;
  using AlphaPoint = Kernel::Point_d;
  using Alpha_complex = Gudhi::alpha_complex::Alpha_complex<Kernel>;

private:
  vector<double> voids_diameter;
  vector<std::pair<Filtration_value, Filtration_value>> voids_persistence;

public:
  /**
   * @brief Default constructor
   *
   */
  AlphaFiltration() = default;

  /**
   * @brief Construct a new Alpha Filtration object using a vector of points
   *
   * This constructor converts the input vector into a vector of AlphaPoint
   * needed for the Gudhi Alpha_complex construction.
   *
   * @param points
   */
  AlphaFiltration(const vector<Point<CoordType>> &points);

  /**
   * @brief Computes an estimate of the maximum diameter of each void using
   * voids persistence information
   *
   */
  void compute_voids_diameter();

  /**
   * @brief Prints the persistence information and th estimate of the maximum
   * diameter of the nmax most persistent voids
   *
   * @param nmax Number od diameters to print
   */
  void print_diameters(unsigned nmax) const;

  /**
   * @brief Calls Filtration::compute_persistent_cohomology and
   * compute_voids_diameter
   *
   */
  void make_analysis() override;

  /**
   * @brief Prints a report of the alpha filtration analysis
   *
   */
  void print_analysis() override;
};

// ##########################################################################

/**
 * @brief Derived class from Filtration. It represents a radial filtration
 *
 */
class RadialFiltration : public Filtration {
public:
  unsigned n_edges;
  double max_radius;
  /// Number of connected components in dimension 0 barcodes with persistence
  /// less or equal the 10% of the max_radius divided by the number of vessel
  /// segments
  double tortuosity_descriptor;
  /// Number of loops divided by the number of vessel segments
  double loops_descriptor;

  /**
   * @brief Default constructor
   *
   */
  RadialFiltration() = default;

  /**
   * @brief Construct a new Radial Filtration object using a vector of points
   * and a vector of edges
   *
   * This contructor determines the centre of mass of the nodes and then assign
   * to each node a filtration value equal to the distance from the centre. Edge
   * filtration values are assigned as the maximum filtration value of the two
   * nodes.
   *
   * @param points
   * @param edges
   */
  RadialFiltration(const vector<Point<CoordType>> &points,
                   const vector<Edge> &edges);

  /**
   * @brief Computes the tortuosity and loops descriptors
   */
  void compute_descriptors();

  /**
   * @brief Calls Filtration::compute_persistent_cohomology and
   * compute_descriptors
   *
   */
  void make_analysis() override;

  /**
   * @brief Prints a report of the radial filtration analysis
   *
   */
  void print_analysis() override;
};

// ##########################################################################

/**
 * @brief Class that parses and stores vertices and edges of a network.
 *
 */
class NetworkData {
public:
  vector<Point<CoordType>> vertices;
  vector<Edge> edges;

  /**
   * @brief Prints the number of vertices and edges of the network.
   *
   */
  void print_network_info() const;

  /**
   * @brief Reads data from files.
   * @param filename_v Name of a file that contains the vertex data for the
   * network (including path)
   * @param sep_v Separator character used in the vertex data file (like ','
   * or ' ')
   * @param filename_e Name of a file that contains the edge data for the
   * network. (including path)
   * @param sep_e Separator character used in the edge data file (like ','
   * or ' ')
   * @param from_0 Boolean flag that indicates whether the vertices indices in
   * the edge data file are 0-based (i.e., start from 0) or 1-based (i.e., start
   * from 1)
   */
  void parse(const std::string &filename_v, char sep_v,
             const string &filename_e, char sep_e, bool from_0);

private:
  /**
   * @brief Reads vertices data from a file
   * @param filename Name of a file that contains the vertex data for the
   * network (including path)
   * @param sep Separator character used in the vertex data file (like ',' or '
   * ')
   */
  void read_vertices(const string &filename, char sep);

  /**
   * @brief Reads edge data from a file
   * @param filename Name of a file that contains the edge data for the
   * network. (including path)
   * @param sep Separator character used in the edge data file (like ',' or ' ')
   * @param from_0 Boolean flag that indicates whether the vertex indices in
   * the edge data file are 0-based (i.e., start from 0) or 1-based (i.e., start
   * from 1)
   */
  void read_edges(const string &filename, char sep, bool from_0);
};

// ##########################################################################

/**
 * @brief Class that performs the topological analysis of a network.
 *
 */
class NetworkAnalysis {
private:
  const NetworkData *data;
  std::vector<std::unique_ptr<Filtration>> filtrations;

public:
  /**
   * @brief Default constructor
   */
  NetworkAnalysis() = default;

  /**
   * @brief Construct a new Network Analysis object building a radial and an
   * alpha filtration.
   *
   * @param data_
   */
  NetworkAnalysis(const NetworkData *data_);

  /**
   * @brief Performs the topological analysis of the network using all the
   * filtrations.
   *
   */
  void analyse() const;

  /**
   * @brief Prints the results of the topological analysis.
   *
   */
  void print_global_analysis() const;

  /**
   * @brief Saves the persistence information of all the filtrations to files.
   *
   */
  void save_persistence(const std::string &name) const;
};