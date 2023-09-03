/**
 * @file geometry.h
 * @author Davide Gurrieri (davide.gurrieri@mail.polimi.it)
 * @brief This file contains the geometry concepts needed to represent a
 * vascular network.
 * @date 2023-08-29
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <vector>
/// Represents an edge of the network
using Edge = std::array<unsigned, 2>;

/**
 * @brief Template class that represents a point in a generic metric space
 *
 * @tparam T Type of the coordinates of the point
 */
template <typename T> class Point {
public:
  std::vector<T> coordinates;

  /**
   * @brief Construct a Point object from a vector of coordinates
   *
   * @param coords
   */
  Point(const std::vector<T> &coords) : coordinates(coords) {}

  /**
   * @brief Construct a default initialized Point object of dimension n
   *
   * @param n
   */
  Point(unsigned int n) : coordinates(std::vector<T>(n, T())) {}

  /**
   * @brief Returns the dimension of the point
   *
   * @return unsigned
   */
  unsigned dimension() const { return coordinates.size(); }

  /**
   * @brief Sum operator
   *
   * @param other
   * @return Point<T>
   */
  Point<T> operator+(const Point<T> &other) const {
    if (coordinates.size() != other.coordinates.size()) {
      throw std::runtime_error("Mismatched dimensions");
    }
    std::vector<T> result;
    result.reserve(coordinates.size());
    for (size_t i = 0; i < coordinates.size(); ++i) {
      result.push_back(coordinates[i] + other.coordinates[i]);
    }
    return Point<T>(result);
  }

  /**
   * @brief Subtraction operator
   *
   * @param other
   * @return Point<T>
   */
  Point<T> operator-(const Point<T> &other) const {
    if (coordinates.size() != other.coordinates.size()) {
      throw std::runtime_error("Mismatched dimensions");
    }
    std::vector<T> result;
    result.reserve(coordinates.size());
    for (size_t i = 0; i < coordinates.size(); ++i) {
      result.push_back(coordinates[i] - other.coordinates[i]);
    }
    return Point<T>(result);
  }

  /**
   * @brief Multiplication by a scalar operator
   *
   * @param scalar
   * @return Point<T>
   */
  Point<T> operator*(const T scalar) const {
    std::vector<T> result;
    result.reserve(coordinates.size());
    for (size_t i = 0; i < coordinates.size(); ++i) {
      result.push_back(coordinates[i] * scalar);
    }
    return Point<T>(result);
  }

  /**
   * @brief Division by a scalar operator
   *
   * @param scalar
   * @return Point<T>
   */
  Point<T> operator/(const T scalar) const {
    if (scalar == 0) {
      throw std::runtime_error("Division by zero");
    }
    std::vector<T> result;
    result.reserve(coordinates.size());
    for (size_t i = 0; i < coordinates.size(); ++i) {
      result.push_back(coordinates[i] / scalar);
    }
    return Point<T>(result);
  }

  /**
   * @brief Returns the Euclidean norm of the point
   *
   * @return T
   */
  T norm() const {
    T sum_of_squares = 0;
    for (const T &coord : coordinates) {
      sum_of_squares += coord * coord;
    }
    return std::sqrt(sum_of_squares);
  }

  /**
   * @brief << operator that prints the coordinates of the point
   *
   * @return std::ostream&
   */
  friend std::ostream &operator<<(std::ostream &os, const Point<T> &point) {
    os << "(";
    for (size_t i = 0; i < point.coordinates.size(); ++i) {
      os << point.coordinates[i];
      if (i < point.coordinates.size() - 1) {
        os << ", ";
      }
    }
    os << ")";
    return os;
  }
};