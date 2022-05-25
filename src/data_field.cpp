// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "src/data_field.hpp"
#include <fftw3.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
// Class type declarations
template class DataField<int>;
template class DataField<double>;
template class DataField<fftw_complex>;
// _____________________________________________________________________________
template <typename T>
void DataField<T>::create(size_t ac, size_t as) {
  this->array_count = ac;
  this->array_size = as;
  for (size_t i = 0; i < array_count; ++i) {
    arrays.push_back(new T[array_size]);
  }
  allocated_memory = true;
}
template <>
void DataField<double>::create(size_t ac, size_t as) {
  this->array_count = ac;
  this->array_size = as;
  for (size_t i = 0; i < array_count; ++i) {
    arrays.push_back(fftw_alloc_real(array_size));
  }
  allocated_memory = true;
}
template <>
void DataField<fftw_complex>::create(size_t ac, size_t as) {
  this->array_count = ac;
  this->array_size = as;
  for (size_t i = 0; i < array_count; ++i) {
    arrays.push_back(fftw_alloc_complex(array_count));
  }
  allocated_memory = true;
}
// _____________________________________________________________________________
template <typename T>
DataField<T>::DataField()
  : array_count(0),
    array_size(0),
    bin_width(1.),
    allocated_memory(false) {
}
// _____________________________________________________________________________
template <typename T>
DataField<T>::DataField(size_t array_count, size_t array_size)
  : array_count(array_count),
    array_size(array_size),
    bin_width(1.),
    allocated_memory(false) {
  this->create(array_count, array_size);
  this->zeros();
}
// _____________________________________________________________________________
template <typename T>
DataField<T>::DataField(size_t array_count, size_t array_size, double bin_width)
  : array_count(array_count),
    array_size(array_size),
    bin_width(bin_width),
    allocated_memory(false) {
  this->create(array_count, array_size);
  this->zeros();
}
// _____________________________________________________________________________
template <typename T>
DataField<T>::DataField(const DataField<T>& other)
  : array_count(other.get_array_count()),
    array_size(other.get_array_size()),
    bin_width(other.get_bin_width()),
    allocated_memory(false) {
  this->create(array_count, array_size);
  // TODO(Moritz): *this = other
  for (size_t i = 0; i < array_count; ++i) {
    for (size_t j = 0; j < array_size; ++j) {
      this->at(i, j) = other.element(i, j);
    }
  }
}
template <>
DataField<fftw_complex>::DataField(
    const DataField<fftw_complex>& other)
  : array_count(other.get_array_count()),
    array_size(other.get_array_size()),
    bin_width(other.get_bin_width()),
    allocated_memory(false) {
  this->create(array_count, array_size);
  // TODO(Moritz): *this = other
  for (size_t i = 0; i < array_count; ++i) {
    for (size_t j = 0; j < array_size; ++j) {
      this->at(i, j)[0] = other.element(i, j)[0];
      this->at(i, j)[1] = other.element(i, j)[1];
    }
  }
}
// _____________________________________________________________________________
template <typename T>
DataField<T>::~DataField() {
  this->clear();
}
// _____________________________________________________________________________
template <typename T>
void DataField<T>::clear() {
  for (auto it = arrays.begin(); it != arrays.end(); ++it) {
    delete [] *it;
  }
  arrays.clear();
  array_count = 0;
  array_size = 0;
  allocated_memory = false;
}
template <>
void DataField<double>::clear() {
  for (auto it = arrays.begin(); it != arrays.end(); ++it) {
    fftw_free(*it);
  }
  arrays.clear();
  array_count = 0;
  array_size = 0;
  allocated_memory = false;
}
template <>
void DataField<fftw_complex>::clear() {
  for (auto it = arrays.begin(); it != arrays.end(); ++it) {
    fftw_free(*it);
  }
  arrays.clear();
  array_count = 0;
  array_size = 0;
  allocated_memory = false;
}
// _____________________________________________________________________________
template <typename T>
void DataField<T>::recreate(size_t ss, size_t as) {
  this->clear();
  this->create(ss, as);
}
// _____________________________________________________________________________
template <typename T>
T& DataField<T>::at(size_t i, size_t j) {
  if (j >= array_size) {
    std::cerr << "DataField::at(): \"Error: Index out of range\".";
    std::cerr << std::endl;
    exit(1);
  }
  return arrays.at(i)[j];
}
// _____________________________________________________________________________
template <typename T>
T& DataField<T>::element(size_t i, size_t j) const {
  return arrays.at(i)[j];
}
// _____________________________________________________________________________
template <typename T>
T* DataField<T>::array(size_t i) {
  return arrays.at(i);
}
// _____________________________________________________________________________
template <typename T>
size_t DataField<T>::get_array_count() const {
  return array_count;
}
// _____________________________________________________________________________
template <typename T>
size_t DataField<T>::get_array_size() const {
  return array_size;
}
// _____________________________________________________________________________
template <typename T>
double DataField<T>::get_bin_width() const {
  return this->bin_width;
}
// _____________________________________________________________________________
template <typename T>
void DataField<T>::set_bin_width(double bw) {
  this->bin_width = bw;
}
// _____________________________________________________________________________
template <typename T>
void DataField<T>::print(
    std::ostream& outstream, std::streamsize stream_size) const {
  // Save previous output precision
  std::streamsize old_stream_size = outstream.precision();
  // Header
  outstream.precision(stream_size);
  outstream << "# r ";
  for (size_t i = 0; i < array_count; ++i) {
    outstream << "A_[" << i << "]";
    outstream << " ";
  }
  outstream << std::endl;
  // Data
  for (size_t i = 0; i < array_size; ++i) {
    outstream << bin_width * static_cast<double>(i+1) << " ";
    for (size_t j = 0; j < array_count; ++j) {
      outstream << this->element(j, i);
      outstream << " ";
    }
    outstream << std::endl;
  }
  // Restore output precision
  outstream.precision(old_stream_size);
}
template <>
void DataField<fftw_complex>::print(
    std::ostream& outstream, std::streamsize stream_size) const {
  // Save previous output precision
  std::streamsize old_stream_size = outstream.precision();
  // Header
  outstream.precision(stream_size);
  outstream << "# r ";
  for (size_t i = 0; i < array_count; ++i) {
    outstream << "A_[" << i << "]";
    outstream << " ";
  }
  outstream << std::endl;
  // Data
  for (size_t i = 0; i < array_size; ++i) {
    outstream << bin_width * static_cast<double>(i+1) << " ";
    for (size_t j = 0; j < array_count; ++j) {
      outstream << this->element(j, i)[0] << "+ i ";  // real part
      outstream << this->element(j, i)[1];  // imaginary part
      outstream << " ";
    }
    outstream << std::endl;
  }
  // Restore output precision
  outstream.precision(old_stream_size);
}
// _____________________________________________________________________________
template <typename T>
void DataField<T>::print(std::ostream& outstream) const {
  print(outstream, std::numeric_limits<double>::max_digits10);
}
// _____________________________________________________________________________
template <typename T>
void DataField<T>::print() const {
  print(std::cout);
}
// _____________________________________________________________________________
template <typename T>
void DataField<T>::set_all_elements_to(const double number) {
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      this->at(i, j) = number;
    }
  }
}
template <>
void DataField<fftw_complex>::set_all_elements_to(
    const double number) {
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      this->at(i, j)[0] = number;
      this->at(i, j)[1] = 0.;
    }
  }
}
// _____________________________________________________________________________
template <typename T>
void DataField<T>::zeros() {
  set_all_elements_to(0.);
}
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
template <typename T>
DataField<T>& DataField<T>::operator=(DataField<T> other) {
  if (this->array_count == other.get_array_count() &&
      this->array_size == other.get_array_size()
  ) {
    for (size_t i = 0; i < this->array_count; ++i) {
      for (size_t j = 0; j < this->array_size; ++j) {
        this->at(i, j) = other.element(i, j);
      }
    }
    this->bin_width = other.get_bin_width();
  } else {
    std::cerr << "DataField::operator=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot copy.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DataField<fftw_complex>& DataField<fftw_complex>::operator=(
    DataField<fftw_complex> other) {
  if (this->array_count == other.get_array_count() &&
      this->array_size == other.get_array_size()
  ) {
    for (size_t i = 0; i < this->array_count; ++i) {
      for (size_t j = 0; j < this->array_size; ++j) {
        this->at(i, j)[0] = other.element(i, j)[0];
        this->at(i, j)[1] = other.element(i, j)[1];
      }
    }
    this->bin_width = other.get_bin_width();
  } else {
    std::cerr << "DataField::operator=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot copy.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
// _____________________________________________________________________________
template <typename T>
DataField<T>& DataField<T>::operator+=(DataField<T> other) {
  if (this->array_count == other.get_array_count() &&
      this->array_size == other.get_array_size()
  ) {
    for (size_t i = 0; i < this->array_count; ++i) {
      for (size_t j = 0; j < this->array_size; ++j) {
        this->at(i, j) += other.element(i, j);
      }
    }
  } else {
    std::cerr << "DataField::operator+=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot add.\"";
    std::cerr << std::endl;
    exit(1);
  }
  if (fabs(this->bin_width - other.get_bin_width()) > 1e-12) {
    std::cerr << "DataField::operator+=():";
    std::cerr << "\"Warning: Matrices do not have the same bin width.\"";
    std::cerr << std::endl;
  }
  return *this;
}
template <>
DataField<fftw_complex>& DataField<fftw_complex>::operator+=(
    DataField<fftw_complex> other) {
  if (this->array_count == other.get_array_count() &&
      this->array_size == other.get_array_size()
  ) {
    for (size_t i = 0; i < this->array_count; ++i) {
      for (size_t j = 0; j < this->array_size; ++j) {
        this->at(i, j)[0] += other.element(i, j)[0];
        this->at(i, j)[1] += other.element(i, j)[1];
      }
    }
  } else {
    std::cerr << "DataField::operator+=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot add.\"";
    std::cerr << std::endl;
    exit(1);
  }
  if (fabs(this->bin_width - other.get_bin_width()) > 1e-12) {
    std::cerr << "DataField::operator+=():";
    std::cerr << "\"Warning: Matrices do not have the same bin width.\"";
    std::cerr << std::endl;
  }
  return *this;
}
// _____________________________________________________________________________
template <typename T>
DataField<T>& DataField<T>::operator-=(DataField<T> other) {
  if (this->array_count == other.get_array_count() &&
      this->array_size == other.get_array_size()
  ) {
    for (size_t i = 0; i < this->array_count; ++i) {
      for (size_t j = 0; j < this->array_size; ++j) {
        this->at(i, j) -= other.element(i, j);
      }
    }
  } else {
    std::cerr << "DataField::operator-=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  if (fabs(this->bin_width - other.get_bin_width()) > 1e-12) {
    std::cerr << "DataField::operator-=():";
    std::cerr << "\"Warning: Matrices do not have the same bin width.\"";
    std::cerr << std::endl;
  }
  return *this;
}
template <>
DataField<fftw_complex>& DataField<fftw_complex>::operator-=(
    DataField<fftw_complex> other) {
  if (this->array_count == other.get_array_count() &&
      this->array_size == other.get_array_size()
  ) {
    for (size_t i = 0; i < this->array_count; ++i) {
      for (size_t j = 0; j < this->array_size; ++j) {
        this->at(i, j)[0] -= other.element(i, j)[0];  // real
        this->at(i, j)[1] -= other.element(i, j)[1];  // imaginary
      }
    }
  } else {
    std::cerr << "DataField::operator-=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  if (fabs(this->bin_width - other.get_bin_width()) > 1e-12) {
    std::cerr << "DataField::operator-=():";
    std::cerr << "\"Warning: Matrices do not have the same bin width.\"";
    std::cerr << std::endl;
  }
  return *this;
}
// _____________________________________________________________________________
template <typename T>
DataField<T>& DataField<T>::operator*=(DataField<T> other) {
  if (this->array_count == other.get_array_count() &&
      this->array_size == other.get_array_size()
  ) {
    for (size_t i = 0; i < this->array_count; ++i) {
      for (size_t j = 0; j < this->array_size; ++j) {
        this->at(i, j) *= other.element(i, j);
      }
    }
  } else {
    std::cerr << "DataField::operator*=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  if (fabs(this->bin_width - other.get_bin_width()) > 1e-12) {
    std::cerr << "DataField::operator*=():";
    std::cerr << "\"Warning: Matrices do not have the same bin width.\"";
    std::cerr << std::endl;
  }
  return *this;
}
template <>
DataField<fftw_complex>& DataField<fftw_complex>::operator*=(
    DataField<fftw_complex> other) {
  if (this->array_count == other.get_array_count() &&
      this->array_size == other.get_array_size()
  ) {
    for (size_t i = 0; i < this->array_count; ++i) {
      for (size_t j = 0; j < this->array_size; ++j) {
        this->at(i, j)[0] *= other.element(i, j)[0];
        this->at(i, j)[1] *= other.element(i, j)[1];
      }
    }
  } else {
    std::cerr << "DataField::operator*=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  if (fabs(this->bin_width - other.get_bin_width()) > 1e-12) {
    std::cerr << "DataField::operator*=():";
    std::cerr << "\"Warning: Matrices do not have the same bin width.\"";
    std::cerr << std::endl;
  }
  return *this;
}
// _____________________________________________________________________________
template <typename T>
DataField<T>& DataField<T>::operator/=(DataField<T> other) {
  if (this->array_count == other.get_array_count() &&
      this->array_size == other.get_array_size()
  ) {
    for (size_t i = 0; i < this->array_count; ++i) {
      for (size_t j = 0; j < this->array_size; ++j) {
        this->at(i, j) /= other.element(i, j);
      }
    }
  } else {
    std::cerr << "DataField::operator/=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  if (fabs(this->bin_width - other.get_bin_width()) > 1e-12) {
    std::cerr << "DataField::operator/=():";
    std::cerr << "\"Warning: Matrices do not have the same bin width.\"";
    std::cerr << std::endl;
  }
  return *this;
}
template <>
DataField<fftw_complex>& DataField<fftw_complex>::operator/=(
    DataField<fftw_complex> other) {
  if (this->array_count == other.get_array_count() &&
      this->array_size == other.get_array_size()
  ) {
    for (size_t i = 0; i < this->array_count; ++i) {
      for (size_t j = 0; j < this->array_size; ++j) {
        this->at(i, j)[0] /= other.element(i, j)[0];
        this->at(i, j)[1] /= other.element(i, j)[1];
      }
    }
  } else {
    std::cerr << "DataField::operator/=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  if (fabs(this->bin_width - other.get_bin_width()) > 1e-12) {
    std::cerr << "DataField::operator/=():";
    std::cerr << "\"Warning: Matrices do not have the same bin width.\"";
    std::cerr << std::endl;
  }
  return *this;
}
// _____________________________________________________________________________
template <typename T>
DataField<T>& DataField<T>::operator+=(double other) {
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      this->at(i, j) += other;
    }
  }
  return *this;
}
template <>
DataField<fftw_complex>& DataField<fftw_complex>::operator+=(
    double other) {
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      this->at(i, j)[0] += other;
    }
  }
  return *this;
}
// _____________________________________________________________________________
template <typename T>
DataField<T>& DataField<T>::operator-=(double other) {
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      this->at(i, j) -= other;
    }
  }
  return *this;
}
template <>
DataField<fftw_complex>& DataField<fftw_complex>::operator-=(
    double other) {
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      this->at(i, j)[0] -= other;
    }
  }
  return *this;
}
// _____________________________________________________________________________
template <typename T>
DataField<T>& DataField<T>::operator*=(double other) {
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      this->at(i, j) *= other;
    }
  }
  return *this;
}
template <>
DataField<fftw_complex>& DataField<fftw_complex>::operator*=(
    double other) {
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      this->at(i, j)[0] *= other;
      this->at(i, j)[1] *= other;
    }
  }
  return *this;
}
// _____________________________________________________________________________
template <typename T>
DataField<T>& DataField<T>::operator/=(double other) {
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      this->at(i, j) /= other;
    }
  }
  return *this;
}
template <>
DataField<fftw_complex>& DataField<fftw_complex>::operator/=(
    double other) {
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      this->at(i, j)[0] /= other;
      this->at(i, j)[1] /= other;
    }
  }
  return *this;
}
// _____________________________________________________________________________
template <typename T>
DataField<T> DataField<T>::operator+(DataField<T> other) {
  DataField<T> result(other.get_array_count(), other.get_array_size());
  if (this->array_count == other.get_array_count() &&
      this->array_size == other.get_array_size()
  ) {
    for (size_t i = 0; i < this->array_count; ++i) {
      for (size_t j = 0; j < this->array_size; ++j) {
        result.at(i, j) = this->element(i, j) + other.element(i, j);
      }
    }
  } else {
    std::cerr << "DataField::operator+():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot add.\"";
    std::cerr << std::endl;
    exit(1);
  }
  if (fabs(this->bin_width - other.get_bin_width()) > 1e-12) {
    std::cerr << "DataField::operator+():";
    std::cerr << "\"Warning: Matrices do not have the same bin width.\"";
    std::cerr << std::endl;
  }
  return result;
}
template <>
DataField<fftw_complex> DataField<fftw_complex>::operator+(
    DataField<fftw_complex> other) {
  DataField<fftw_complex> result(
      other.get_array_count(), other.get_array_size());
  if (this->array_count == other.get_array_count() &&
      this->array_size == other.get_array_size()
  ) {
    for (size_t i = 0; i < this->array_count; ++i) {
      for (size_t j = 0; j < this->array_size; ++j) {
        result.at(i, j)[0] = this->element(i, j)[0] + other.element(i, j)[0];
        result.at(i, j)[1] = this->element(i, j)[1] + other.element(i, j)[1];
      }
    }
  } else {
    std::cerr << "DataField::operator+():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot add.\"";
    std::cerr << std::endl;
    exit(1);
  }
  if (fabs(this->bin_width - other.get_bin_width()) > 1e-12) {
    std::cerr << "DataField::operator+():";
    std::cerr << "\"Warning: Matrices do not have the same bin width.\"";
    std::cerr << std::endl;
  }
  return result;
}
// _____________________________________________________________________________
template <typename T>
DataField<T> DataField<T>::operator-(DataField<T> other) {
  DataField<T> result(other.get_array_count(), other.get_array_size());
  if (this->array_count == other.get_array_count() &&
      this->array_size == other.get_array_size()
  ) {
    for (size_t i = 0; i < this->array_count; ++i) {
      for (size_t j = 0; j < this->array_size; ++j) {
        result.at(i, j) = this->element(i, j) - other.element(i, j);
      }
    }
  } else {
    std::cerr << "DataField::operator-():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  if (fabs(this->bin_width - other.get_bin_width()) > 1e-12) {
    std::cerr << "DataField::operator-():";
    std::cerr << "\"Warning: Matrices do not have the same bin width.\"";
    std::cerr << std::endl;
  }
  return result;
}
template <>
DataField<fftw_complex> DataField<fftw_complex>::operator-(
    DataField<fftw_complex> other) {
  DataField<fftw_complex> result(
      other.get_array_count(), other.get_array_size());
  if (this->array_count == other.get_array_count() &&
      this->array_size == other.get_array_size()
  ) {
    for (size_t i = 0; i < this->array_count; ++i) {
      for (size_t j = 0; j < this->array_size; ++j) {
        result.at(i, j)[0] = this->element(i, j)[0] - other.element(i, j)[0];
        result.at(i, j)[1] = this->element(i, j)[1] - other.element(i, j)[1];
      }
    }
  } else {
    std::cerr << "DataField::operator-():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  if (fabs(this->bin_width - other.get_bin_width()) > 1e-12) {
    std::cerr << "DataField::operator-():";
    std::cerr << "\"Warning: Matrices do not have the same bin width.\"";
    std::cerr << std::endl;
  }
  return result;
}
// _____________________________________________________________________________
template <typename T>
DataField<T> DataField<T>::operator*(DataField<T> other) {
  DataField<T> result(other.get_array_count(), other.get_array_size());
  if (this->array_count == other.get_array_count() &&
      this->array_size == other.get_array_size()
  ) {
    for (size_t i = 0; i < this->array_count; ++i) {
      for (size_t j = 0; j < this->array_size; ++j) {
        result.at(i, j) = this->element(i, j) * other.element(i, j);
      }
    }
  } else {
    std::cerr << "DataField::operator*():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot add.\"";
    std::cerr << std::endl;
    exit(1);
  }
  if (fabs(this->bin_width - other.get_bin_width()) > 1e-12) {
    std::cerr << "DataField::operator*():";
    std::cerr << "\"Warning: Matrices do not have the same bin width.\"";
    std::cerr << std::endl;
  }
  return result;
}
template <>
DataField<fftw_complex> DataField<fftw_complex>::operator*(
    DataField<fftw_complex> other) {
  DataField<fftw_complex> result(
      other.get_array_count(), other.get_array_size());
  if (this->array_count == other.get_array_count() &&
      this->array_size == other.get_array_size()
  ) {
    for (size_t i = 0; i < this->array_count; ++i) {
      for (size_t j = 0; j < this->array_size; ++j) {
        result.at(i, j)[0] = this->element(i, j)[0] * other.element(i, j)[0] -
            this->element(i, j)[1] * other.element(i, j)[1];
        result.at(i, j)[1] = this->element(i, j)[0] * other.element(i, j)[1] +
            this->element(i, j)[1] * other.element(i, j)[0];
      }
    }
  } else {
    std::cerr << "DataField::operator*():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot add.\"";
    std::cerr << std::endl;
    exit(1);
  }
  if (fabs(this->bin_width - other.get_bin_width()) > 1e-12) {
    std::cerr << "DataField::operator*():";
    std::cerr << "\"Warning: Matrices do not have the same bin width.\"";
    std::cerr << std::endl;
  }
  return result;
}
// _____________________________________________________________________________
template <typename T>
DataField<T> DataField<T>::operator/(DataField<T> other) {
  DataField<T> result(other.get_array_count(), other.get_array_size());
  if (this->array_count == other.get_array_count() &&
      this->array_size == other.get_array_size()
  ) {
    for (size_t i = 0; i < this->array_count; ++i) {
      for (size_t j = 0; j < this->array_size; ++j) {
        result.at(i, j) = this->element(i, j) / other.element(i, j);
      }
    }
  } else {
    std::cerr << "DataField::operator/():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  if (fabs(this->bin_width - other.get_bin_width()) > 1e-12) {
    std::cerr << "DataField::operator/():";
    std::cerr << "\"Warning: Matrices do not have the same bin width.\"";
    std::cerr << std::endl;
  }
  return result;
}
template <>
DataField<fftw_complex> DataField<fftw_complex>::operator/(
    DataField<fftw_complex> other) {
  DataField<fftw_complex> result(
      other.get_array_count(), other.get_array_size());
  if (this->array_count == other.get_array_count() &&
      this->array_size == other.get_array_size()
  ) {
    for (size_t i = 0; i < this->array_count; ++i) {
      for (size_t j = 0; j < this->array_size; ++j) {
        result.at(i, j)[0] = (this->element(i, j)[0] * other.element(i, j)[0] +
            this->element(i, j)[1] * other.element(i, j)[1]) /
            (other.element(i, j)[0] * other.element(i, j)[0] +
            other.element(i, j)[1] * other.element(i, j)[1]);
        result.at(i, j)[1] = (this->element(i, j)[1] * other.element(i, j)[0] -
            this->element(i, j)[0] * other.element(i, j)[1]) /
            (other.element(i, j)[0] * other.element(i, j)[0] +
            other.element(i, j)[1] * other.element(i, j)[1]);
      }
    }
  } else {
    std::cerr << "DataField::operator/():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  if (fabs(this->bin_width - other.get_bin_width()) > 1e-12) {
    std::cerr << "DataField::operator/():";
    std::cerr << "\"Warning: Matrices do not have the same bin width.\"";
    std::cerr << std::endl;
  }
  return result;
}
// _____________________________________________________________________________
template <typename T>
DataField<T> DataField<T>::operator+(double other) {
  DataField<T> result(this->array_count, this->array_size);
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      result.at(i, j) = this->element(i, j) + other;
    }
  }
  return result;
}
template <>
DataField<fftw_complex> DataField<fftw_complex>::operator+(
    double other) {
  DataField<fftw_complex> result(this->array_count, this->array_size);
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      result.at(i, j)[0] = this->element(i, j)[0] + other;
      result.at(i, j)[1] = this->element(i, j)[1];
    }
  }
  return result;
}
// _____________________________________________________________________________
template <typename T>
DataField<T> DataField<T>::operator-(double other) {
  DataField<T> result(this->array_count, this->array_size);
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      result.at(i, j) = this->element(i, j) - other;
    }
  }
  return result;
}
template <>
DataField<fftw_complex> DataField<fftw_complex>::operator-(
    double other) {
  DataField<fftw_complex> result(this->array_count, this->array_size);
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      result.at(i, j)[0] = this->element(i, j)[0] - other;
      result.at(i, j)[1] = this->element(i, j)[1];
    }
  }
  return result;
}
// _____________________________________________________________________________
template <typename T>
DataField<T> DataField<T>::operator*(double other) {
  DataField<T> result(this->array_count, this->array_size);
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      result.at(i, j) = this->element(i, j) * other;
    }
  }
  return result;
}
template <>
DataField<fftw_complex> DataField<fftw_complex>::operator*(
    double other) {
  DataField<fftw_complex> result(this->array_count, this->array_size);
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      result.at(i, j)[0] = this->element(i, j)[0] * other;
      result.at(i, j)[1] = this->element(i, j)[1] * other;
    }
  }
  return result;
}
// _____________________________________________________________________________
template <typename T>
DataField<T> DataField<T>::operator/(double other) {
  DataField<T> result(this->array_count, this->array_size);
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      result.at(i, j) = this->element(i, j) / other;
    }
  }
  return result;
}
template <>
DataField<fftw_complex> DataField<fftw_complex>::operator/(
    double other) {
  DataField<fftw_complex> result(this->array_count, this->array_size);
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      result.at(i, j)[0] = this->element(i, j)[0] / other;
      result.at(i, j)[1] = this->element(i, j)[1] / other;
    }
  }
  return result;
}
// _____________________________________________________________________________
template <typename U>
DataField<U> operator+(double current, DataField<U> other) {  // friend
  DataField<U> result(other.get_array_count(), other.get_array_size());
  for (size_t i = 0; i < other.get_array_count(); ++i) {
    for (size_t j = 0; j < other.get_array_size(); ++j) {
      result.at(i, j) = current + other.element(i, j);
    }
  }
  return result;
}
template <>
DataField<fftw_complex> operator+(
    double current, DataField<fftw_complex> other) {  // friend
  DataField<fftw_complex> result(
      other.get_array_count(), other.get_array_size());
  for (size_t i = 0; i < other.get_array_count(); ++i) {
    for (size_t j = 0; j < other.get_array_size(); ++j) {
      result.at(i, j)[0] = current + other.element(i, j)[0];
      result.at(i, j)[1] = other.element(i, j)[1];
    }
  }
  return result;
}
// _____________________________________________________________________________
template <typename V>
DataField<V> operator-(double current, DataField<V> other) {  // friend
  DataField<V> result(other.get_array_count(), other.get_array_size());
  for (size_t i = 0; i < other.get_array_count(); ++i) {
    for (size_t j = 0; j < other.get_array_size(); ++j) {
      result.at(i, j) = current - other.element(i, j);
    }
  }
  return result;
}
template <>
DataField<fftw_complex> operator-(
    double current, DataField<fftw_complex> other) {  // friend
  DataField<fftw_complex> result(
      other.get_array_count(), other.get_array_size());
  for (size_t i = 0; i < other.get_array_count(); ++i) {
    for (size_t j = 0; j < other.get_array_size(); ++j) {
      result.at(i, j)[0] = current - other.element(i, j)[0];
      result.at(i, j)[1] = - other.element(i, j)[1];
    }
  }
  return result;
}
// _____________________________________________________________________________
template <typename W>
DataField<W> operator*(double current, DataField<W> other) {  // friend
  DataField<W> result(other.get_array_count(), other.get_array_size());
  for (size_t i = 0; i < other.get_array_count(); ++i) {
    for (size_t j = 0; j < other.get_array_size(); ++j) {
      result.at(i, j) = current * other.element(i, j);
    }
  }
  return result;
}
template <>
DataField<fftw_complex> operator*(
    double current, DataField<fftw_complex> other) {  // friend
  DataField<fftw_complex> result(
      other.get_array_count(), other.get_array_size());
  for (size_t i = 0; i < other.get_array_count(); ++i) {
    for (size_t j = 0; j < other.get_array_size(); ++j) {
      result.at(i, j)[0] = current * other.element(i, j)[0];
      result.at(i, j)[1] = current * other.element(i, j)[1];
    }
  }
  return result;
}
// _____________________________________________________________________________
template <typename X>
DataField<X> operator/(double current, DataField<X> other) {  // friend
  DataField<X> result(other.get_array_count(), other.get_array_size());
  for (size_t i = 0; i < other.get_array_count(); ++i) {
    for (size_t j = 0; j < other.get_array_size(); ++j) {
      result.at(i, j) = current / other.element(i, j);
    }
  }
  return result;
}
template <>
DataField<fftw_complex> operator/(
    double current, DataField<fftw_complex> other) {  // friend
  DataField<fftw_complex> result(
      other.get_array_count(), other.get_array_size());
  for (size_t i = 0; i < other.get_array_count(); ++i) {
    for (size_t j = 0; j < other.get_array_size(); ++j) {
      result.at(i, j)[0] = current * other.element(i, j)[0] /
          (other.element(i, j)[0] * other.element(i, j)[0] +
          other.element(i, j)[1] * other.element(i, j)[1]);
      result.at(i, j)[1] = -current * other.element(i, j)[1] /
          (other.element(i, j)[0] * other.element(i, j)[0] +
          other.element(i, j)[1] * other.element(i, j)[1]);
    }
  }
  return result;
}
// _____________________________________________________________________________
