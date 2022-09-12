// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2019 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "data_frame.hpp"  // NOLINT
#include <fftw3.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
// Template class forward declarations
template class DataFrame<1, double>;
template class DataFrame<2, double>;
template class DataFrame<3, double>;
template class DataFrame<1, fftw_complex>;
template class DataFrame<2, fftw_complex>;
template class DataFrame<3, fftw_complex>;
// _____________________________________________________________________________
template <size_t dim, typename T>
size_t DataFrame<dim, T>::size() const {
  return array_size;
}
// _____________________________________________________________________________
template <size_t dim, typename T>
std::vector<size_t> DataFrame<dim, T>::size_dim() const {
  return array_dimensions;
}
// _____________________________________________________________________________
template <size_t dim, typename T>
T& DataFrame<dim, T>::at(size_t i) {
  if (i >= array_size) {
    std::cerr << "DataFrame::at(): \"Error: Index out of range.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return data[i];
}
// _____________________________________________________________________________
template <size_t dim, typename T>
T& DataFrame<dim, T>::at(size_t i, size_t j) {
  std::cerr << "DataFrame::at(): \"Error: This accessor only works in 2D.\"";
  std::cerr << std::endl;
  exit(1);
  return data[0];
}
template <>
double& DataFrame<2, double>::at(size_t i, size_t j) {
  return this->at(coordinates_to_index(std::vector<size_t>{i, j}));
}
template <>
fftw_complex& DataFrame<2, fftw_complex>::at(size_t i, size_t j) {
  return this->at(coordinates_to_index(std::vector<size_t>{i, j}));
}
// _____________________________________________________________________________
template <size_t dim, typename T>
T& DataFrame<dim, T>::at(size_t i, size_t j, size_t k) {
  std::cerr << "DataFrame::at(): \"Error: This accessor only works in 3D.\"";
  std::cerr << std::endl;
  exit(1);
  return data[0];
}
template <>
double& DataFrame<3, double>::at(size_t i, size_t j, size_t k) {
  return this->at(coordinates_to_index(std::vector<size_t>{i, j, k}));
}
template <>
fftw_complex& DataFrame<3, fftw_complex>::at(
    size_t i, size_t j, size_t k) {
  return this->at(coordinates_to_index(std::vector<size_t>{i, j, k}));
}
// _____________________________________________________________________________
template <size_t dim, typename T>
T& DataFrame<dim, T>::element(size_t i) const {
  if (i >= array_size) {
    std::cerr << "DataFrame::element():";
    std::cerr << "\"Error: Index out of range.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return data[i];
}
// _____________________________________________________________________________
template <size_t dim, typename T>
T& DataFrame<dim, T>::element(size_t i, size_t j) const {
  std::cerr << "DataFrame::element(): ";
  std::cerr << "\"Error: This accessor only works in 2D.\"";
  std::cerr << std::endl;
  exit(1);
  return data[0];
}
template <>
double& DataFrame<2, double>::element(size_t i, size_t j) const {
  return this->element(coordinates_to_index(std::vector<size_t>{i, j}));
}
template <>
fftw_complex& DataFrame<2, fftw_complex>::element(size_t i, size_t j) const {
  return this->element(coordinates_to_index(std::vector<size_t>{i, j}));
}
// _____________________________________________________________________________
template <size_t dim, typename T>
T& DataFrame<dim, T>::element(size_t i, size_t j, size_t k) const {
  std::cerr << "DataFrame::element(): ";
  std::cerr << "\"Error: This accessor only works in 3D.\"";
  std::cerr << std::endl;
  exit(1);
  return data[0];
}
template <>
double& DataFrame<3, double>::element(size_t i, size_t j, size_t k) const {
  return this->element(coordinates_to_index(std::vector<size_t>{i, j, k}));
}
template <>
fftw_complex& DataFrame<3, fftw_complex>::element(
    size_t i, size_t j, size_t k) const {
  return this->element(coordinates_to_index(std::vector<size_t>{i, j, k}));
}
// _____________________________________________________________________________
template <size_t dim, typename T>
std::string DataFrame<dim, T>::element_string(
    size_t i, std::streamsize stream_size) const {
  if (i >= array_size) {
    std::cerr << "DataFrame::element():";
    std::cerr << "\"Error: Index out of range.\"";
    std::cerr << std::endl;
    exit(1);
  }
  std::stringstream string_stream;
  string_stream.precision(stream_size);
  string_stream << data[i];
  return string_stream.str();
}
template <>
std::string DataFrame<1, fftw_complex>::element_string(
    size_t i, std::streamsize stream_size) const {
  if (i >= array_size) {
    std::cerr << "DataFrame::element():";
    std::cerr << "\"Error: Index out of range.\"";
    std::cerr << std::endl;
    exit(1);
  }
  std::ostringstream string_stream;
  string_stream.precision(stream_size);
  string_stream << data[i][0] << "+i" << data[i][1];
  return string_stream.str();
}
template <>
std::string DataFrame<2, fftw_complex>::element_string(
    size_t i, std::streamsize stream_size) const {
  if (i >= array_size) {
    std::cerr << "DataFrame::element():";
    std::cerr << "\"Error: Index out of range.\"";
    std::cerr << std::endl;
    exit(1);
  }
  std::ostringstream string_stream;
  string_stream.precision(stream_size);
  string_stream << data[i][0] << "+i" << data[i][1];
  return string_stream.str();
}
template <>
std::string DataFrame<3, fftw_complex>::element_string(
    size_t i, std::streamsize stream_size) const {
  if (i >= array_size) {
    std::cerr << "DataFrame::element():";
    std::cerr << "\"Error: Index out of range.\"";
    std::cerr << std::endl;
    exit(1);
  }
  std::ostringstream string_stream;
  string_stream.precision(stream_size);
  string_stream << data[i][0] << "+i" << data[i][1];
  return string_stream.str();
}
// _____________________________________________________________________________
template <size_t dim, typename T>
T* DataFrame<dim, T>::array() {
  return data;
}
// _____________________________________________________________________________
template <size_t dim, typename T>
void DataFrame<dim, T>::set_all_elements_to(T value) {
  for (size_t i = 0; i < array_size; ++i) {
    this->at(i) = value;
  }
}
template <>
void DataFrame<1, fftw_complex>::set_all_elements_to(fftw_complex value) {
  for (size_t i = 0; i < array_size; ++i) {
    this->at(i)[0] = value[0];
    this->at(i)[1] = value[1];
  }
}
template <>
void DataFrame<2, fftw_complex>::set_all_elements_to(fftw_complex value) {
  for (size_t i = 0; i < array_size; ++i) {
    this->at(i)[0] = value[0];
    this->at(i)[1] = value[1];
  }
}
template <>
void DataFrame<3, fftw_complex>::set_all_elements_to(fftw_complex value) {
  for (size_t i = 0; i < array_size; ++i) {
    this->at(i)[0] = value[0];
    this->at(i)[1] = value[1];
  }
}
// _____________________________________________________________________________
template <size_t dim, typename T>
void DataFrame<dim, T>::zero() {
  set_all_elements_to(0.);
}
template <>
void DataFrame<1, fftw_complex>::zero() {
  fftw_complex zero{0., 0.};
  set_all_elements_to(zero);
}
template <>
void DataFrame<2, fftw_complex>::zero() {
  fftw_complex zero{0., 0.};
  set_all_elements_to(zero);
}
template <>
void DataFrame<3, fftw_complex>::zero() {
  fftw_complex zero{0., 0.};
  set_all_elements_to(zero);
}
// _____________________________________________________________________________
template <size_t dim, typename T>
DataFrame<dim, T>::DataFrame()
  : array_size(0),
    array_dimensions(std::vector<size_t>(dim, 0)) {
  data = new T[array_size];
}
// _____________________________________________________________________________
template <size_t dim, typename T>
DataFrame<dim, T>::DataFrame(std::vector<size_t> array_dimensions)
  : array_size(0),
    array_dimensions(array_dimensions) {
  if (array_dimensions.size() != dim) {
    std::cerr << "DataFrame::DataFrame(): ";
    std::cerr << "\"Error: Dimensions of DataFrame and constructor";
    std::cerr << "argument do not match.\"" << std::endl;
    exit(1);
  }
  array_size = 1;
  for (auto it = array_dimensions.begin(); it != array_dimensions.end(); ++it) {
    array_size *= *it;
  }
  data = new T[array_size];
  zero();
}
// _____________________________________________________________________________
template <size_t dim, typename T>
DataFrame<dim, T>::DataFrame(size_t array_size) : DataFrame() {
  std::cerr << "DataFrame::DataFrame(): ";
  std::cerr << "\"Error: Specifying just the array_size only works in one-";
  std::cerr << "dimensional DataFrames.\"" << std::endl;
  exit(1);
}
template <>
DataFrame<1, double>::DataFrame(size_t array_size)
  : DataFrame<1, double>(std::vector<size_t>(1, array_size)) {
}
template <>
DataFrame<1, fftw_complex>::DataFrame(size_t array_size)
  : DataFrame<1, fftw_complex>(std::vector<size_t>(1, array_size)) {
}
// _____________________________________________________________________________
template <size_t dim, typename T>
DataFrame<dim, T>::DataFrame(const DataFrame<dim, T>& other)
  : array_size(other.size()),
    array_dimensions(other.size_dim()) {
  data = new T[array_size];
  for (size_t i = 0; i < array_size; ++i) {
    this->at(i) = other.element(i);
  }
}
template <>
DataFrame<1, fftw_complex>::DataFrame(
    const DataFrame<1, fftw_complex>& other)
  : array_size(other.size()),
    array_dimensions(other.size_dim()) {
  data = new fftw_complex[array_size];
  for (size_t i = 0; i < array_size; ++i) {
    this->at(i)[0] = other.element(i)[0];  // Real
    this->at(i)[1] = other.element(i)[1];  // Imaginary
  }
}
template <>
DataFrame<2, fftw_complex>::DataFrame(
    const DataFrame<2, fftw_complex>& other)
  : array_size(other.size()),
    array_dimensions(other.size_dim()) {
  data = new fftw_complex[array_size];
  for (size_t i = 0; i < array_size; ++i) {
    this->at(i)[0] = other.element(i)[0];  // Real
    this->at(i)[1] = other.element(i)[1];  // Imaginary
  }
}
template <>
DataFrame<3, fftw_complex>::DataFrame(
    const DataFrame<3, fftw_complex>& other)
  : array_size(other.size()),
    array_dimensions(other.size_dim()) {
  data = new fftw_complex[array_size];
  for (size_t i = 0; i < array_size; ++i) {
    this->at(i)[0] = other.element(i)[0];  // Real
    this->at(i)[1] = other.element(i)[1];  // Imaginary
  }
}
// _____________________________________________________________________________
template <size_t dim, typename T>
DataFrame<dim, T>::~DataFrame() {
  delete [] data;
}
// _____________________________________________________________________________
template <size_t dim, typename T>
bool DataFrame<dim, T>::same_size(const DataFrame<dim, T>& other) const {
  return (this->size() == other.size());
}
// _____________________________________________________________________________
template <size_t dim, typename T>
void DataFrame<dim, T>::print(
    std::ostream& outstream, std::streamsize stream_size) const {
  // Save previous output precision
  std::streamsize old_stream_size = outstream.precision();
  // Header
  outstream.precision(stream_size);
  outstream << "# [index] [data]" << std::endl;
  // Data
  std::vector<size_t> coordinates;
  for (size_t i = 0; i < array_size; ++i) {
    coordinates = index_to_coordinates(i);
    outstream << "(";
    for (auto it = coordinates.begin(); it != coordinates.end(); ++it) {
      outstream << *it;
      if (it != (coordinates.end() - 1)) {
        outstream << ",";
      } else {
        outstream << ")";
      }
    }
    outstream << " " << this->element_string(i, stream_size) << std::endl;
  }
  // Restore output precision
  outstream.precision(old_stream_size);
}
// _____________________________________________________________________________
template <size_t dim, typename T>
void DataFrame<dim, T>::print(std::ostream& outstream) const {
  print(outstream, std::numeric_limits<double>::max_digits10);
}
// _____________________________________________________________________________
template <size_t dim, typename T>
void DataFrame<dim, T>::print() const {
  print(std::cout);
}
// _____________________________________________________________________________
template <size_t dim, typename T>
DataFrame<dim, T>& DataFrame<dim, T>::operator=(
    const DataFrame<dim, T>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) = other.element(i);
    }
  } else {
    std::cerr << "DataFrame::operator=():";
    std::cerr << " \"ERROR: Array sizes do not match. Cannot copy.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrame<1, fftw_complex>& DataFrame<1, fftw_complex>::operator=(
    const DataFrame<1, fftw_complex>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i)[0] = other.element(i)[0];  // Real
      this->at(i)[1] = other.element(i)[1];  // Imaginary
    }
  } else {
    std::cerr << "DataFrame::operator=():";
    std::cerr << " \"ERROR: Array sizes do not match. Cannot copy.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrame<2, fftw_complex>& DataFrame<2, fftw_complex>::operator=(
    const DataFrame<2, fftw_complex>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i)[0] = other.element(i)[0];  // Real
      this->at(i)[1] = other.element(i)[1];  // Imaginary
    }
  } else {
    std::cerr << "DataFrame::operator=():";
    std::cerr << " \"ERROR: Array sizes do not match. Cannot copy.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrame<3, fftw_complex>& DataFrame<3, fftw_complex>::operator=(
    const DataFrame<3, fftw_complex>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i)[0] = other.element(i)[0];  // Real
      this->at(i)[1] = other.element(i)[1];  // Imaginary
    }
  } else {
    std::cerr << "DataFrame::operator=():";
    std::cerr << " \"ERROR: Array sizes do not match. Cannot copy.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
// _____________________________________________________________________________
template <size_t dim, typename T>
DataFrame<dim, T>& DataFrame<dim, T>::operator+=(
    const DataFrame<dim, T>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) += other.element(i);
    }
  } else {
    std::cerr << "DataFrame::operator+=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot add.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrame<1, fftw_complex>& DataFrame<1, fftw_complex>::operator+=(
    const DataFrame<1, fftw_complex>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i)[0] += other.element(i)[0];  // Real
      this->at(i)[1] += other.element(i)[1];  // Imaginary
    }
  } else {
    std::cerr << "DataFrame::operator+=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot add.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrame<2, fftw_complex>& DataFrame<2, fftw_complex>::operator+=(
    const DataFrame<2, fftw_complex>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i)[0] += other.element(i)[0];  // Real
      this->at(i)[1] += other.element(i)[1];  // Imaginary
    }
  } else {
    std::cerr << "DataFrame::operator+=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot add.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrame<3, fftw_complex>& DataFrame<3, fftw_complex>::operator+=(
    const DataFrame<3, fftw_complex>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i)[0] += other.element(i)[0];  // Real
      this->at(i)[1] += other.element(i)[1];  // Imaginary
    }
  } else {
    std::cerr << "DataFrame::operator+=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot add.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
// _____________________________________________________________________________
template <size_t dim, typename T>
DataFrame<dim, T>& DataFrame<dim, T>::operator-=(
    const DataFrame<dim, T>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) -= other.element(i);
    }
  } else {
    std::cerr << "DataFrame::operator-=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrame<1, fftw_complex>& DataFrame<1, fftw_complex>::operator-=(
    const DataFrame<1, fftw_complex>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i)[0] -= other.element(i)[0];  // Real
      this->at(i)[1] -= other.element(i)[1];  // Imaginary
    }
  } else {
    std::cerr << "DataFrame::operator-=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrame<2, fftw_complex>& DataFrame<2, fftw_complex>::operator-=(
    const DataFrame<2, fftw_complex>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i)[0] -= other.element(i)[0];  // Real
      this->at(i)[1] -= other.element(i)[1];  // Imaginary
    }
  } else {
    std::cerr << "DataFrame::operator-=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrame<3, fftw_complex>& DataFrame<3, fftw_complex>::operator-=(
    const DataFrame<3, fftw_complex>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i)[0] -= other.element(i)[0];  // Real
      this->at(i)[1] -= other.element(i)[1];  // Imaginary
    }
  } else {
    std::cerr << "DataFrame::operator-=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
// _____________________________________________________________________________
template <size_t dim, typename T>
DataFrame<dim, T>& DataFrame<dim, T>::operator*=(
    const DataFrame<dim, T>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) *= other.element(i);
    }
  } else {
    std::cerr << "DataFrame::operator*=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot multiply.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrame<1, fftw_complex>& DataFrame<1, fftw_complex>::operator*=(
    const DataFrame<1, fftw_complex>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i)[0] = this->element(i)[0] * other.element(i)[0] -
          this->element(i)[1] * other.element(i)[1];  // Real
      this->at(i)[1] = this->element(i)[0] * other.element(i)[1] +
          this->element(i)[1] * other.element(i)[0];  // Imaginary
    }
  } else {
    std::cerr << "DataFrame::operator*=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot multiply.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrame<2, fftw_complex>& DataFrame<2, fftw_complex>::operator*=(
    const DataFrame<2, fftw_complex>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i)[0] = this->element(i)[0] * other.element(i)[0] -
          this->element(i)[1] * other.element(i)[1];  // Real
      this->at(i)[1] = this->element(i)[0] * other.element(i)[1] +
          this->element(i)[1] * other.element(i)[0];  // Imaginary
    }
  } else {
    std::cerr << "DataFrame::operator*=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot multiply.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrame<3, fftw_complex>& DataFrame<3, fftw_complex>::operator*=(
    const DataFrame<3, fftw_complex>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i)[0] = this->element(i)[0] * other.element(i)[0] -
          this->element(i)[1] * other.element(i)[1];  // Real
      this->at(i)[1] = this->element(i)[0] * other.element(i)[1] +
          this->element(i)[1] * other.element(i)[0];  // Imaginary
    }
  } else {
    std::cerr << "DataFrame::operator*=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot multiply.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
// _____________________________________________________________________________
template <size_t dim, typename T>
DataFrame<dim, T>& DataFrame<dim, T>::operator/=(
    const DataFrame<dim, T>& other) {
  if (this->size() == other.size()) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) /= other.element(i);
    }
  } else {
    std::cerr << "DataFrame::operator/=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot divide.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrame<1, fftw_complex>& DataFrame<1, fftw_complex>::operator/=(
    const DataFrame<1, fftw_complex>& other) {
  if (this->size() == other.size()) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i)[0] = (this->element(i)[0] * other.element(i)[0] +
          this->element(i)[1] * other.element(i)[1]) /
          (pow(other.element(i)[0], 2) + pow(other.element(i)[1], 2));  // Real
      this->at(i)[1] = (this->element(i)[1] * other.element(i)[0] -
          this->element(i)[0] * other.element(i)[1]) /
          (pow(other.element(i)[0], 2) + pow(other.element(i)[1], 2));  // Imag
    }
  } else {
    std::cerr << "DataFrame::operator/=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot divide.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrame<2, fftw_complex>& DataFrame<2, fftw_complex>::operator/=(
    const DataFrame<2, fftw_complex>& other) {
  if (this->size() == other.size()) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i)[0] = (this->element(i)[0] * other.element(i)[0] +
          this->element(i)[1] * other.element(i)[1]) /
          (pow(other.element(i)[0], 2) + pow(other.element(i)[1], 2));  // Real
      this->at(i)[1] = (this->element(i)[1] * other.element(i)[0] -
          this->element(i)[0] * other.element(i)[1]) /
          (pow(other.element(i)[0], 2) + pow(other.element(i)[1], 2));  // Imag
    }
  } else {
    std::cerr << "DataFrame::operator/=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot divide.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrame<3, fftw_complex>& DataFrame<3, fftw_complex>::operator/=(
    const DataFrame<3, fftw_complex>& other) {
  if (this->size() == other.size()) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i)[0] = (this->element(i)[0] * other.element(i)[0] +
          this->element(i)[1] * other.element(i)[1]) /
          (pow(other.element(i)[0], 2) + pow(other.element(i)[1], 2));  // Real
      this->at(i)[1] = (this->element(i)[1] * other.element(i)[0] -
          this->element(i)[0] * other.element(i)[1]) /
          (pow(other.element(i)[0], 2) + pow(other.element(i)[1], 2));  // Imag
    }
  } else {
    std::cerr << "DataFrame::operator/=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot divide.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
// _____________________________________________________________________________
template <size_t dim, typename T>
DataFrame<dim, T>& DataFrame<dim, T>::operator*=(const double other) {
  for (size_t i = 0; i < this->array_size; ++i) {
    this->at(i) *= other;
  }
  return *this;
}
template <>
DataFrame<1, fftw_complex>& DataFrame<1, fftw_complex>::operator*=(
    const double other) {
  for (size_t i = 0; i < this->array_size; ++i) {
    this->at(i)[0] = this->element(i)[0] * other;  // Real
    this->at(i)[1] = this->element(i)[1] * other;  // Imaginary
  }
  return *this;
}
template <>
DataFrame<2, fftw_complex>& DataFrame<2, fftw_complex>::operator*=(
    const double other) {
  for (size_t i = 0; i < this->array_size; ++i) {
    this->at(i)[0] = this->element(i)[0] * other;  // Real
    this->at(i)[1] = this->element(i)[1] * other;  // Imaginary
  }
  return *this;
}
template <>
DataFrame<3, fftw_complex>& DataFrame<3, fftw_complex>::operator*=(
    const double other) {
  for (size_t i = 0; i < this->array_size; ++i) {
    this->at(i)[0] = this->element(i)[0] * other;  // Real
    this->at(i)[1] = this->element(i)[1] * other;  // Imaginary
  }
  return *this;
}
// _____________________________________________________________________________
template <size_t dim, typename T>
DataFrame<dim, T> DataFrame<dim, T>::operator+(
    const DataFrame<dim, T>& other) {
  DataFrame<dim, T> result(*this);
  if (this->same_size(other)) {
    result += other;
  } else {
    std::cerr << "DataFrame::operator+():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot add.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return result;
}
// _____________________________________________________________________________
template <size_t dim, typename T>
DataFrame<dim, T> DataFrame<dim, T>::operator-(
    const DataFrame<dim, T>& other) {
  DataFrame<dim, T> result(*this);
  if (this->same_size(other)) {
    result -= other;
  } else {
    std::cerr << "DataFrame::operator-():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return result;
}
// _____________________________________________________________________________
template <size_t dim, typename T>
DataFrame<dim, T> DataFrame<dim, T>::operator*(
    const DataFrame<dim, T>& other) {
  DataFrame<dim, T> result(*this);
  if (this->same_size(other)) {
    result *= other;
  } else {
    std::cerr << "DataFrame::operator*():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot multiply.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return result;
}
// _____________________________________________________________________________
template <size_t dim, typename T>
DataFrame<dim, T> DataFrame<dim, T>::operator/(
    const DataFrame<dim, T>& other) {
  DataFrame<dim, T> result(*this);
  if (this->same_size(other)) {
    result /= other;
  } else {
    std::cerr << "DataFrame::operator/():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot divide.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return result;
}
// _____________________________________________________________________________
template <size_t dim, typename T>
DataFrame<dim, T> DataFrame<dim, T>::operator*(const double other) {
  DataFrame<dim, T> result(*this);
  result *= other;
  return result;
}
// _____________________________________________________________________________
DataFrame<1, double> operator*(const double current,
    const DataFrame<1, double>& other) {  // friend
  DataFrame<1, double> result(other);
  result *= current;
  return result;
}
DataFrame<2, double> operator*(const double current,
    const DataFrame<2, double>& other) {  // friend
  DataFrame<2, double> result(other);
  result *= current;
  return result;
}
DataFrame<3, double> operator*(const double current,
    const DataFrame<3, double>& other) {  // friend
  DataFrame<3, double> result(other);
  result *= current;
  return result;
}
DataFrame<1, fftw_complex> operator*(const double current,
    const DataFrame<1, fftw_complex>& other) {  // friend
  DataFrame<1, fftw_complex> result(other);
  result *= current;
  return result;
}
DataFrame<2, fftw_complex> operator*(const double current,
    const DataFrame<2, fftw_complex>& other) {  // friend
  DataFrame<2, fftw_complex> result(other);
  result *= current;
  return result;
}
DataFrame<3, fftw_complex> operator*(const double current,
    const DataFrame<3, fftw_complex>& other) {  // friend
  DataFrame<3, fftw_complex> result(other);
  result *= current;
  return result;
}
// _____________________________________________________________________________
DataFrame<1, double> exp(const DataFrame<1, double>& other) {  // friend
  DataFrame<1, double> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i) = exp(other.element(i));
  }
  return result;
}
DataFrame<2, double> exp(const DataFrame<2, double>& other) {  // friend
  DataFrame<2, double> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i) = exp(other.element(i));
  }
  return result;
}
DataFrame<3, double> exp(const DataFrame<3, double>& other) {  // friend
  DataFrame<3, double> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i) = exp(other.element(i));
  }
  return result;
}
DataFrame<1, fftw_complex> exp(const DataFrame<1, fftw_complex>& other) {  // friend  // NOLINT
  DataFrame<1, fftw_complex> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i)[0] = exp(other.element(i)[0]) * cos(other.element(i)[1]);
    result.at(i)[1] = exp(other.element(i)[0]) * sin(other.element(i)[1]);
  }
  return result;
}
DataFrame<2, fftw_complex> exp(const DataFrame<2, fftw_complex>& other) {  // friend  // NOLINT
  DataFrame<2, fftw_complex> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i)[0] = exp(other.element(i)[0]) * cos(other.element(i)[1]);
    result.at(i)[1] = exp(other.element(i)[0]) * sin(other.element(i)[1]);
  }
  return result;
}
DataFrame<3, fftw_complex> exp(const DataFrame<3, fftw_complex>& other) {  // friend  // NOLINT
  DataFrame<3, fftw_complex> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i)[0] = exp(other.element(i)[0]) * cos(other.element(i)[1]);
    result.at(i)[1] = exp(other.element(i)[0]) * sin(other.element(i)[1]);
  }
  return result;
}
// _____________________________________________________________________________
DataFrame<1, double> log_natural(const DataFrame<1, double>& other) {  // friend
  DataFrame<1, double> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i) = log(other.element(i));
  }
  return result;
}
DataFrame<2, double> log_natural(const DataFrame<2, double>& other) {  // friend
  DataFrame<2, double> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i) = log(other.element(i));
  }
  return result;
}
DataFrame<3, double> log_natural(const DataFrame<3, double>& other) {  // friend
  DataFrame<3, double> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i) = log(other.element(i));
  }
  return result;
}
DataFrame<1, fftw_complex> log_natural(
    const DataFrame<1, fftw_complex>& other) {  // friend
  double re{0.}, im{0.}, radius{0.}, angle{0.};
  DataFrame<1, fftw_complex> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    re = other.element(i)[0];
    im = other.element(i)[1];
    radius = sqrt(re * re + im * im);
    if (radius == 0.)  {
      angle = 0.;
    } else if (re > 0. && im >= 0.) {
      angle = asin(im / radius);
    } else if (re <= 0. && im >= 0.) {
      angle = .5 * M_PI + asin(fabs(re) / radius);
    } else if (re <= 0. && im < 0.) {
      angle = M_PI + asin(fabs(im) / radius);
    } else if (re > 0. && im < 0.) {
      angle = 1.5 * M_PI + asin(re / radius);
    }
    result.at(i)[0] = log(radius);
    result.at(i)[1] = angle;
  }
  std::cerr << "DataFrame::log_natural():";
  std::cerr << " \"Warning: The complex logarithm is used.";
  std::cerr << "Keep in mind, that there are infinite solutions.\"";
  std::cerr << std::endl;
  return result;
}
DataFrame<2, fftw_complex> log_natural(
    const DataFrame<2, fftw_complex>& other) {  // friend
  DataFrame<1, fftw_complex> other_1d(other.size());
  for (size_t i = 0; i < other.size(); ++i) {
    other_1d.at(i)[0] = other.element(i)[0];
    other_1d.at(i)[1] = other.element(i)[1];
  }
  other_1d = log_natural(other_1d);
  DataFrame<2, fftw_complex> result(other);
  for (size_t i = 0; i < other.size(); ++i) {
    result.at(i)[0] = other_1d.at(i)[0];
    result.at(i)[1] = other_1d.at(i)[1];
  }
  return result;
}
DataFrame<3, fftw_complex> log_natural(
    const DataFrame<3, fftw_complex>& other) {  // friend
  DataFrame<1, fftw_complex> other_1d(other.size());
  for (size_t i = 0; i < other.size(); ++i) {
    other_1d.at(i)[0] = other.element(i)[0];
    other_1d.at(i)[1] = other.element(i)[1];
  }
  other_1d = log_natural(other_1d);
  DataFrame<3, fftw_complex> result(other);
  for (size_t i = 0; i < other.size(); ++i) {
    result.at(i)[0] = other_1d.at(i)[0];
    result.at(i)[1] = other_1d.at(i)[1];
  }
  return result;
}
// _____________________________________________________________________________
DataFrame<1, double> abs(const DataFrame<1, double>& other) {  // friend
  DataFrame<1, double> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i) = fabs(result.at(i));
  }
  return result;
}
DataFrame<2, double> abs(const DataFrame<2, double>& other) {  // friend
  DataFrame<2, double> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i) = fabs(result.at(i));
  }
  return result;
}
DataFrame<3, double> abs(const DataFrame<3, double>& other) {  // friend
  DataFrame<3, double> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i) = fabs(result.at(i));
  }
  return result;
}
DataFrame<1, fftw_complex> abs(const DataFrame<1, fftw_complex>& other) {  // friend  // NOLINT
  DataFrame<1, fftw_complex> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i)[0] = sqrt(pow(result.at(i)[0], 2) + pow(result.at(i)[1], 2));
    result.at(i)[1] = 0.;
  }
  return result;
}
DataFrame<2, fftw_complex> abs(const DataFrame<2, fftw_complex>& other) {  // friend  // NOLINT
  DataFrame<2, fftw_complex> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i)[0] = sqrt(pow(result.at(i)[0], 2) + pow(result.at(i)[1], 2));
    result.at(i)[1] = 0.;
  }
  return result;
}
DataFrame<3, fftw_complex> abs(const DataFrame<3, fftw_complex>& other) {  // friend  // NOLINT
  DataFrame<3, fftw_complex> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i)[0] = sqrt(pow(result.at(i)[0], 2) + pow(result.at(i)[1], 2));
    result.at(i)[1] = 0.;
  }
  return result;
}
// _____________________________________________________________________________
double max(const DataFrame<1, double>& other) {  // friend
  return *std::max_element(other.data, other.data + other.array_size);
}
double max(const DataFrame<2, double>& other) {  // friend
  return *std::max_element(other.data, other.data + other.array_size);
}
double max(const DataFrame<3, double>& other) {  // friend
  return *std::max_element(other.data, other.data + other.array_size);
}
// _____________________________________________________________________________
template <size_t dim, typename T>
std::vector<size_t> DataFrame<dim, T>::index_to_coordinates(size_t index)
    const {
  std::vector<size_t> coordinates{};
  std::vector<size_t> array_dims = array_dimensions;
  array_dims.resize(3, 1);  // project dim!=3 into dim==3
  coordinates.push_back(
      (index % (array_dims.at(0) * array_dims.at(1))) % array_dims.at(0));
  coordinates.push_back(
      (index % (array_dims.at(0) * array_dims.at(1))) / array_dims.at(0));
  coordinates.push_back(index / (array_dims.at(0) * array_dims.at(1)));
  coordinates.resize(dim);  // remove trivial coordinates from projection
  return coordinates;
}
// _____________________________________________________________________________
template <size_t dim, typename T>
size_t DataFrame<dim, T>::coordinates_to_index(
    std::vector<size_t> coordinates) const {
  if (array_dimensions.size() != coordinates.size()) {
    std::cerr << "DataFrame::coordinates_to_index(): ";
    std::cerr << "\"Error: Wrong number of dimensions supplied.\"";
    std::cerr << std::endl;
    exit(1);
  }
  size_t index{0};
  std::vector<size_t> bins{1};  // bins that need to be added for diff. dims
  for (size_t i = 0; i < coordinates.size(); ++i) {
    index += coordinates.at(i) * bins.at(i);
    bins.push_back(bins.at(i) * array_dimensions.at(i));
  }
  return index;
}
// _____________________________________________________________________________
