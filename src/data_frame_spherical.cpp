// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2022 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "src/data_frame_spherical.hpp"
#include "src/data_frame.hpp"
#include <fftw3.h>
// Class type declarations
template class DataFrameSpherical<double>;
// template class DataFrameSpherical<int>;
// template class DataFrameSpherical<fftw_complex>;
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>::DataFrameSpherical(const Properties properties)
  : DataFrame() {
  // TODO(Moritz): properties --> memory allocation, initialization, ...
}
template <>
DataFrameSpherical<double>::DataFrameSpherical(const Properties properties)
  : DataFrame() {
  // TODO(Moritz): properties --> memory allocation, initialization, ...
}
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>::DataFrameSpherical(const DataFrameSpherical<T>& other) {
  // TODO(Moritz): implement copy constructor
}
template <>
DataFrameSpherical<double>::DataFrameSpherical(
    const DataFrameSpherical<double>& other) {
  // TODO(Moritz): implement copy constructor
}
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>::~DataFrameSpherical() {
}
template <>
DataFrameSpherical<double>::~DataFrameSpherical() {
}
// _____________________________________________________________________________
template <typename T>
size_t DataFrameSpherical<T>::size() const {
  return array_size;
}
// _____________________________________________________________________________
template <typename T>
bool DataFrameSpherical<T>::same_size(const DataFrameSpherical<T>& other)
    const {
  return (this->size() == other.size());
}
// _____________________________________________________________________________
template <typename T>
T& DataFrameSpherical<T>::at(size_t i) {
  if (i >= array_size) {
    std::cerr << "DataFrameSpherical::at(): \"Error: Index out of range\".";
    std::cerr << std::endl;
    exit(1);
  }
  return data[i];
}
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>& DataFrameSpherical<T>::operator=(const DataFrameSpherical<T>& other) {
  if (this->size() == other.size()) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) = other.at(i);
    }
  } else {
    std::cerr << "DataFrameSpherical::operator=():";
    std::cerr << " \"ERROR: Array sizes do not match. Cannot copy.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
// TODO(Moritz): implement a constant version of at() like element()
// TODO(Moritz): initialize memory for array in constructor
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>& DataFrameSpherical<T>::operator+=(const DataFrameSpherical<T>& other) {
  //
}
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>& DataFrameSpherical<T>::operator-=(const DataFrameSpherical<T>& other) {
  //
}
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>& DataFrameSpherical<T>::operator*=(const DataFrameSpherical<T>& other) {
  //
}
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>& DataFrameSpherical<T>::operator/=(const DataFrameSpherical<T>& other) {
  //
}
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>& DataFrameSpherical<T>::operator*=(const double other) {
  //
}
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
