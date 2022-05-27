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
size_t DataFrameSpherical<T>::size() const {
  return array_size;
}
template <>
size_t DataFrameSpherical<double>::size() const {
  return array_size;
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
template <>
double& DataFrameSpherical<double>::at(size_t i) {
  if (i >= array_size) {
    std::cerr << "DataFrameSpherical::at(): \"Error: Index out of range\".";
    std::cerr << std::endl;
    exit(1);
  }
  return data[i];
}
// _____________________________________________________________________________
template <typename T>
T& DataFrameSpherical<T>::element(size_t i) const {
  if (i >= array_size) {
    std::cerr << "DataFrameSpherical::element():";
    std::cerr << "\"Error: Index out of range\".";
    std::cerr << std::endl;
    exit(1);
  }
  return data[i];
}
template <>
double& DataFrameSpherical<double>::element(size_t i) const {
  if (i >= array_size) {
    std::cerr << "DataFrameSpherical::element():";
    std::cerr << "\"Error: Index out of range\".";
    std::cerr << std::endl;
    exit(1);
  }
  return data[i];
}
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>::DataFrameSpherical(const Properties& properties)
  : DataFrame() {
  properties.get_property("grid count", &array_size);
  data = new T[array_size];
}
template <>
DataFrameSpherical<double>::DataFrameSpherical(const Properties& properties)
  : DataFrame() {
  properties.get_property("grid count", &array_size);
  data = new double[array_size];
}
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>::DataFrameSpherical(const DataFrameSpherical<T>& other) {
  this->array_size = other.size();
  data = new T[array_size];
  for (size_t i = 0; i < array_size; ++i) {
    this->at(i) = other.element(i);
  }
}
template <>
DataFrameSpherical<double>::DataFrameSpherical(
    const DataFrameSpherical<double>& other) {
  this->array_size = other.size();
  data = new double[array_size];
  for (size_t i = 0; i < array_size; ++i) {
    this->at(i) = other.element(i);
  }
}
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>::~DataFrameSpherical() {
  delete [] data;
}
template <>
DataFrameSpherical<double>::~DataFrameSpherical() {
  delete [] data;
}
// _____________________________________________________________________________
template <typename T>
bool DataFrameSpherical<T>::same_size(const DataFrameSpherical<T>& other)
    const {
  return (this->size() == other.size());
}
template <>
bool DataFrameSpherical<double>::same_size(
    const DataFrameSpherical<double>& other) const {
  return (this->size() == other.size());
}
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>& DataFrameSpherical<T>::operator=(
    const DataFrameSpherical<T>& other) {
  if (this->size() == other.size()) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) = other.element(i);
    }
  } else {
    std::cerr << "DataFrameSpherical::operator=():";
    std::cerr << " \"ERROR: Array sizes do not match. Cannot copy.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrameSpherical<double>& DataFrameSpherical<double>::operator=(
    const DataFrameSpherical<double>& other) {
  if (this->size() == other.size()) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) = other.element(i);
    }
  } else {
    std::cerr << "DataFrameSpherical::operator=():";
    std::cerr << " \"ERROR: Array sizes do not match. Cannot copy.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>& DataFrameSpherical<T>::operator+=(
    const DataFrameSpherical<T>& other) {
  if (this->size() == other.size()) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) += other.element(i);
    }
  } else {
    std::cerr << "DataFrameSpherical::operator+=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot add.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrameSpherical<double>& DataFrameSpherical<double>::operator+=(
    const DataFrameSpherical<double>& other) {
  if (this->size() == other.size()) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) += other.element(i);
    }
  } else {
    std::cerr << "DataFrameSpherical::operator+=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot add.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>& DataFrameSpherical<T>::operator-=(
    const DataFrameSpherical<T>& other) {
  if (this->size() == other.size()) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) -= other.element(i);
    }
  } else {
    std::cerr << "DataFrameSpherical::operator-=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrameSpherical<double>& DataFrameSpherical<double>::operator-=(
    const DataFrameSpherical<double>& other) {
  if (this->size() == other.size()) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) -= other.element(i);
    }
  } else {
    std::cerr << "DataFrameSpherical::operator-=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>& DataFrameSpherical<T>::operator*=(
    const DataFrameSpherical<T>& other) {
  if (this->size() == other.size()) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) *= other.element(i);
    }
  } else {
    std::cerr << "DataFrameSpherical::operator*=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot multiply.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrameSpherical<double>& DataFrameSpherical<double>::operator*=(
    const DataFrameSpherical<double>& other) {
  if (this->size() == other.size()) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) *= other.element(i);
    }
  } else {
    std::cerr << "DataFrameSpherical::operator*=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot multiply.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>& DataFrameSpherical<T>::operator/=(
    const DataFrameSpherical<T>& other) {
  if (this->size() == other.size()) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) /= other.element(i);
    }
  } else {
    std::cerr << "DataFrameSpherical::operator/=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot divide.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrameSpherical<double>& DataFrameSpherical<double>::operator/=(
    const DataFrameSpherical<double>& other) {
  if (this->size() == other.size()) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) /= other.element(i);
    }
  } else {
    std::cerr << "DataFrameSpherical::operator/=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot divide.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>& DataFrameSpherical<T>::operator*=(const double other) {
  for (size_t i = 0; i < this->array_size; ++i) {
    this->at(i) *= other;
  }
  return *this;
}
template <>
DataFrameSpherical<double>& DataFrameSpherical<double>::operator*=(
    const double other) {
  for (size_t i = 0; i < this->array_size; ++i) {
    this->at(i) *= other;
  }
  return *this;
}
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>& DataFrameSpherical<T>::operator+(
    const DataFrameSpherical<T>& other) {
  DataFrameSpherical<T> result(*this);
  //if (this->size() == other.size()) {
  //  result += other;
  //} else {
  //  std::cerr << "DataField::operator+():";
  //  std::cerr << " \"ERROR: Dimensions do not match. Cannot add.\"";
  //  std::cerr << std::endl;
  //  exit(1);
  //}
  //return result;
}
template <>
DataFrameSpherical<double>& DataFrameSpherical<double>::operator+(
    const DataFrameSpherical<double>& other) {
  DataFrameSpherical<double> result(*this);
  //if (this->size() == other.size()) {
  //  result += other;
  //} else {
  //  std::cerr << "DataField::operator+():";
  //  std::cerr << " \"ERROR: Dimensions do not match. Cannot add.\"";
  //  std::cerr << std::endl;
  //  exit(1);
  //}
  //return result;
}
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>& DataFrameSpherical<T>::operator-(
    const DataFrameSpherical<T>& other) {
//  DataFrameSpherical<T> result(*this);
//  if (this->size() == other.size()) {
//    result -= other;
//  } else {
//    std::cerr << "DataField::operator-():";
//    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
//    std::cerr << std::endl;
//    exit(1);
//  }
//  return result;
}
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>& DataFrameSpherical<T>::operator*(
    const DataFrameSpherical<T>& other) {
//  DataFrameSpherical<T> result(*this);
//  if (this->size() == other.size()) {
//    result *= other;
//  } else {
//    std::cerr << "DataField::operator*():";
//    std::cerr << " \"ERROR: Dimensions do not match. Cannot multiply.\"";
//    std::cerr << std::endl;
//    exit(1);
//  }
//  return result;
}
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>& DataFrameSpherical<T>::operator/(
    const DataFrameSpherical<T>& other) {
//  DataFrameSpherical<T> result(*this);
//  if (this->size() == other.size()) {
//    result /= other;
//  } else {
//    std::cerr << "DataField::operator/():";
//    std::cerr << " \"ERROR: Dimensions do not match. Cannot divide.\"";
//    std::cerr << std::endl;
//    exit(1);
//  }
//  return result;
}
// _____________________________________________________________________________
template <typename T>
DataFrameSpherical<T>& DataFrameSpherical<T>::operator*(const double other) {

}
// _____________________________________________________________________________
DataFrameSpherical<T>& DataFrameSpherical<T>::exp() {

}
// _____________________________________________________________________________
DataFrameSpherical<T>& DataFrameSpherical<T>::log_natural() {

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
