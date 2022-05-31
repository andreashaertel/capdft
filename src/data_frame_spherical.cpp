// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2022 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "src/data_frame_spherical.hpp"
#include "src/data_frame.hpp"
#include <fftw3.h>
#include <cmath>
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
// TODO(Moritz): Check if type is correct before casting
template <typename T>
DataFrameSpherical<T>::DataFrameSpherical(const DataFrame& other) {
  const DataFrameSpherical<T>& other_casted =
      dynamic_cast<const DataFrameSpherical<T>&>(other);
  this->array_size = other_casted.size();
  data = new T[array_size];
  for (size_t i = 0; i < array_size; ++i) {
    this->at(i) = other_casted.element(i);
  }
}
template <>
DataFrameSpherical<double>::DataFrameSpherical(const DataFrame& other) {
  const DataFrameSpherical<double>& other_casted =
      dynamic_cast<const DataFrameSpherical<double>&>(other);
  this->array_size = other_casted.size();
  data = new double[array_size];
  for (size_t i = 0; i < array_size; ++i) {
    this->at(i) = other_casted.element(i);
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
    const DataFrame& other) {
  const DataFrameSpherical<T>& other_casted =
      dynamic_cast<const DataFrameSpherical<T>&>(other);
  if (this->same_size(other_casted)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) = other_casted.element(i);
    }
  } else {
    throw bad_size_error;
    std::cerr << "DataFrameSpherical::operator=():";
    std::cerr << " \"ERROR: Array sizes do not match. Cannot copy.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrameSpherical<double>& DataFrameSpherical<double>::operator=(
    const DataFrame& other) {
  const DataFrameSpherical<double>& other_casted =
      dynamic_cast<const DataFrameSpherical<double>&>(other);
  if (this->same_size(other_casted)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) = other_casted.element(i);
    }
  } else {
    throw bad_size_error;
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
    const DataFrame& other) {
  const DataFrameSpherical<T>& other_casted =
      dynamic_cast<const DataFrameSpherical<T>&>(other);
  if (this->same_size(other_casted)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) += other_casted.element(i);
    }
  } else {
    throw bad_size_error;
    std::cerr << "DataFrameSpherical::operator+=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot add.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrameSpherical<double>& DataFrameSpherical<double>::operator+=(
    const DataFrame& other) {
  const DataFrameSpherical<double>& other_casted =
      dynamic_cast<const DataFrameSpherical<double>&>(other);
  if (this->same_size(other_casted)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) += other_casted.element(i);
    }
  } else {
    throw bad_size_error;
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
    const DataFrame& other) {
  const DataFrameSpherical<T>& other_casted =
      dynamic_cast<const DataFrameSpherical<T>&>(other);
  if (this->same_size(other_casted)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) -= other_casted.element(i);
    }
  } else {
    throw bad_size_error;
    std::cerr << "DataFrameSpherical::operator-=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrameSpherical<double>& DataFrameSpherical<double>::operator-=(
    const DataFrame& other) {
  const DataFrameSpherical<double>& other_casted =
      dynamic_cast<const DataFrameSpherical<double>&>(other);
  if (this->same_size(other_casted)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) -= other_casted.element(i);
    }
  } else {
    throw bad_size_error;
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
    const DataFrame& other) {
  const DataFrameSpherical<T>& other_casted =
      dynamic_cast<const DataFrameSpherical<T>&>(other);
  if (this->same_size(other_casted)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) *= other_casted.element(i);
    }
  } else {
    throw bad_size_error;
    std::cerr << "DataFrameSpherical::operator*=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot multiply.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrameSpherical<double>& DataFrameSpherical<double>::operator*=(
    const DataFrame& other) {
  const DataFrameSpherical<double>& other_casted =
      dynamic_cast<const DataFrameSpherical<double>&>(other);
  if (this->same_size(other_casted)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) *= other_casted.element(i);
    }
  } else {
    throw bad_size_error;
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
    const DataFrame& other) {
  const DataFrameSpherical<T>& other_casted =
      dynamic_cast<const DataFrameSpherical<T>&>(other);
  if (this->size() == other_casted.size()) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) /= other_casted.element(i);
    }
  } else {
    throw bad_size_error;
    std::cerr << "DataFrameSpherical::operator/=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot divide.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataFrameSpherical<double>& DataFrameSpherical<double>::operator/=(
    const DataFrame& other) {
  const DataFrameSpherical<double>& other_casted =
      dynamic_cast<const DataFrameSpherical<double>&>(other);
  if (this->size() == other_casted.size()) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) /= other_casted.element(i);
    }
  } else {
    throw bad_size_error;
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
DataFrame DataFrameSpherical<T>::operator+(
    const DataFrame& other) {
  const DataFrameSpherical<T>& other_casted =
      dynamic_cast<const DataFrameSpherical<T>&>(other);
  DataFrameSpherical<T> result(*this);
  if (this->same_size(other_casted)) {
    result += other_casted;
  } else {
    throw bad_size_error;
    std::cerr << "DataFrameSpherical::operator+():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot add.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return result;
}
template <>
DataFrame DataFrameSpherical<double>::operator+(
    const DataFrame& other) {
  const DataFrameSpherical<double>& other_casted =
      dynamic_cast<const DataFrameSpherical<double>&>(other);
  DataFrameSpherical<double> result(*this);
  if (this->same_size(other_casted)) {
    result += other_casted;
  } else {
    throw bad_size_error;
    std::cerr << "DataFrameSpherical::operator+():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot add.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return result;
}
// _____________________________________________________________________________
template <typename T>
DataFrame DataFrameSpherical<T>::operator-(
    const DataFrame& other) {
  const DataFrameSpherical<T>& other_casted =
      dynamic_cast<const DataFrameSpherical<T>&>(other);
  DataFrameSpherical<T> result(*this);
  if (this->same_size(other_casted)) {
    result -= other_casted;
  } else {
    throw bad_size_error;
    std::cerr << "DataFrameSpherical::operator-():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return result;
}
template <>
DataFrame DataFrameSpherical<double>::operator-(
    const DataFrame& other) {
  const DataFrameSpherical<double>& other_casted =
      dynamic_cast<const DataFrameSpherical<double>&>(other);
  DataFrameSpherical<double> result(*this);
  if (this->same_size(other_casted)) {
    result -= other_casted;
  } else {
    throw bad_size_error;
    std::cerr << "DataFrameSpherical::operator-():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return result;
}
// _____________________________________________________________________________
template <typename T>
DataFrame DataFrameSpherical<T>::operator*(
    const DataFrame& other) {
  const DataFrameSpherical<T>& other_casted =
      dynamic_cast<const DataFrameSpherical<T>&>(other);
  DataFrameSpherical<T> result(*this);
  if (this->same_size(other_casted)) {
    result *= other_casted;
  } else {
    throw bad_size_error;
    std::cerr << "DataFrameSpherical::operator*():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot multiply.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return result;
}
template <>
DataFrame DataFrameSpherical<double>::operator*(
    const DataFrame& other) {
  const DataFrameSpherical<double>& other_casted =
      dynamic_cast<const DataFrameSpherical<double>&>(other);
  DataFrameSpherical<double> result(*this);
  if (this->same_size(other_casted)) {
    result *= other_casted;
  } else {
    throw bad_size_error;
    std::cerr << "DataFrameSpherical::operator*():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot multiply.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return result;
}
// _____________________________________________________________________________
template <typename T>
DataFrame DataFrameSpherical<T>::operator/(
    const DataFrame& other) {
  const DataFrameSpherical<T>& other_casted =
      dynamic_cast<const DataFrameSpherical<T>&>(other);
  DataFrameSpherical<T> result(*this);
  if (this->same_size(other_casted)) {
    result /= other_casted;
  } else {
    throw bad_size_error;
    std::cerr << "DataFrameSpherical::operator/():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot divide.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return result;
}
template <>
DataFrame DataFrameSpherical<double>::operator/(
    const DataFrame& other) {
  const DataFrameSpherical<double>& other_casted =
      dynamic_cast<const DataFrameSpherical<double>&>(other);
  DataFrameSpherical<double> result(*this);
  if (this->same_size(other_casted)) {
    result /= other_casted;
  } else {
    throw bad_size_error;
    std::cerr << "DataFrameSpherical::operator/():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot divide.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return result;
}
// _____________________________________________________________________________
template <typename T>
DataFrame DataFrameSpherical<T>::operator*(const T other) {
  DataFrameSpherical<T> result(*this);
  result *= other;
  return result;
}
template <>
DataFrame DataFrameSpherical<double>::operator*(const double other) {
  DataFrameSpherical<double> result(*this);
  result *= other;
  return result;
}
// _____________________________________________________________________________
template <typename U>
DataFrameSpherical<U> operator*(const U current,
    const DataFrameSpherical<U>& other) {
  DataFrameSpherical<U> result(other);
  result *= current;
  return result;
}
template <>
DataFrameSpherical<double> operator*(const double current,
    const DataFrameSpherical<double>& other) {
  DataFrameSpherical<double> result(other);
  result *= current;
  return result;
}
// _____________________________________________________________________________
template <typename V>
DataFrameSpherical<V> exp(DataFrameSpherical<V>& other) {
  DataFrameSpherical<V> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i) = exp(other.element(i));
  }
  return result;
}
template <>
DataFrameSpherical<double> exp(DataFrameSpherical<double>& other) {
  DataFrameSpherical<double> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i) = exp(other.element(i));
  }
  return result;
}
// _____________________________________________________________________________
template <typename W>
DataFrameSpherical<W> log_natural(DataFrameSpherical<W>& other) {
  DataFrameSpherical<W> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i) = log(other.element(i));
  }
  return result;
}
template <>
DataFrameSpherical<double> log_natural(DataFrameSpherical<double>& other) {
  DataFrameSpherical<double> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i) = log(other.element(i));
  }
  return result;
}
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
