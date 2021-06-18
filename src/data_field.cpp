// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "data_field.hpp"
#include <fftw3.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
// _____________________________________________________________________________
DataField::DataField()
  : array_count(0),
    array_size(0),
    bin_width(1.) {
}
// _____________________________________________________________________________
DataField::DataField(size_t array_count, size_t array_size)
  : array_count(array_count),
    array_size(array_size),
    bin_width(1.) {
  this->create(array_count, array_size);
  this->zeros();
}
// _____________________________________________________________________________
DataField::DataField(size_t array_count, size_t array_size, double bin_width)
  : array_count(array_count),
    array_size(array_size),
    bin_width(bin_width) {
  this->create(array_count, array_size);
  this->zeros();
}
// _____________________________________________________________________________
DataField::DataField(const DataField& other)
  : array_count(other.get_array_count()),
    array_size(other.get_array_size()),
    bin_width(other.get_bin_width()) {
  this->create(array_count, array_size);
  for (size_t i = 0; i < array_count; ++i) {
    for (size_t j = 0; j < array_size; ++j) {
      this->at(i, j) = other.element(i, j);
    }
  }
}
// _____________________________________________________________________________
DataField::~DataField() {
  this->clear();
}
// _____________________________________________________________________________
void DataField::clear() {
  for (auto it = arrays.begin(); it != arrays.end(); ++it) {
    fftw_free(*it);
  }
  arrays.clear();
  array_count = 0;
  array_size = 0;
}
// _____________________________________________________________________________
void DataField::create(size_t ac, size_t as) {
  this->array_count = ac;
  this->array_size = as;
  for (size_t i = 0; i < array_count; ++i) {
    arrays.push_back(
        (double*) fftw_malloc(sizeof(double) * array_size));  // NOLINT
  }
}
// _____________________________________________________________________________
void DataField::recreate(size_t ss, size_t as) {
  this->clear();
  this->create(ss, as);
}
// _____________________________________________________________________________
double& DataField::at(size_t i, size_t j) {
  if (j >= array_size) {
    std::cerr << "DataField::at(): \"Error: Index out of range\".";
    std::cerr << std::endl;
    exit(1);
  }
  return arrays.at(i)[j];
}
// _____________________________________________________________________________
double& DataField::element(size_t i, size_t j) const {
  return arrays.at(i)[j];
}
// _____________________________________________________________________________
double* DataField::array(size_t i) {
  return arrays.at(i);
}
// _____________________________________________________________________________
std::vector<double*>::iterator DataField::begin() {
  return arrays.begin();
}
// _____________________________________________________________________________
std::vector<double*>::iterator DataField::end() {
  return arrays.end();
}
// _____________________________________________________________________________
size_t DataField::get_array_count() const {
  return array_count;
}
// _____________________________________________________________________________
size_t DataField::get_array_size() const {
  return array_size;
}
// _____________________________________________________________________________
double DataField::get_bin_width() const {
  return this->bin_width;
}
// _____________________________________________________________________________
void DataField::set_bin_width(double bw) {
  this->bin_width = bw;
}
// _____________________________________________________________________________
void DataField::print(std::ostream& outstream, std::streamsize stream_size)
    const {
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
// _____________________________________________________________________________
void DataField::print(std::ostream& outstream) const {
  print(outstream, std::numeric_limits<double>::max_digits10);
}
// _____________________________________________________________________________
void DataField::print() const {
  print(std::cout);
}
// _____________________________________________________________________________
void DataField::set_all_elements_to(const double number) {
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      this->at(i, j) = number;
    }
  }
}
// _____________________________________________________________________________
void DataField::zeros() {
  set_all_elements_to(0.);
}
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
DataField& DataField::operator=(DataField other) {
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
// _____________________________________________________________________________
DataField& DataField::operator*=(double other) {
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      this->at(i, j) *= other;
    }
  }
  return *this;
}
// _____________________________________________________________________________
DataField& DataField::operator/=(double other) {
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      this->at(i, j) /= other;
    }
  }
  return *this;
}
// _____________________________________________________________________________
DataField& DataField::operator+=(DataField other) {
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
// _____________________________________________________________________________
DataField& DataField::operator-=(DataField other) {
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
// _____________________________________________________________________________
DataField DataField::operator+(DataField other) {
  DataField result(other.get_array_count(), other.get_array_size());
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
// _____________________________________________________________________________
DataField DataField::operator-(DataField other) {
  DataField result(other.get_array_count(), other.get_array_size());
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
// _____________________________________________________________________________
DataField DataField::operator+(double other) {
  DataField result(this->array_count, this->array_size);
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      result.at(i, j) = this->element(i, j) + other;
    }
  }
  return result;
}
// _____________________________________________________________________________
DataField DataField::operator-(double other) {
  DataField result(this->array_count, this->array_size);
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      result.at(i, j) = this->element(i, j) - other;
    }
  }
  return result;
}
// _____________________________________________________________________________
DataField DataField::operator*(double other) {
  DataField result(this->array_count, this->array_size);
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      result.at(i, j) = this->element(i, j) * other;
    }
  }
  return result;
}
// _____________________________________________________________________________
DataField DataField::operator/(double other) {
  DataField result(this->array_count, this->array_size);
  for (size_t i = 0; i < this->array_count; ++i) {
    for (size_t j = 0; j < this->array_size; ++j) {
      result.at(i, j) = this->element(i, j) / other;
    }
  }
  return result;
}
// _____________________________________________________________________________
DataField operator+(double current, DataField other) {  // friend
  DataField result(other.get_array_count(), other.get_array_size());
  for (size_t i = 0; i < other.get_array_count(); ++i) {
    for (size_t j = 0; j < other.get_array_size(); ++j) {
      result.at(i, j) = current + other.element(i, j);
    }
  }
  return result;
}
// _____________________________________________________________________________
DataField operator-(double current, DataField other) {  // friend
  DataField result(other.get_array_count(), other.get_array_size());
  for (size_t i = 0; i < other.get_array_count(); ++i) {
    for (size_t j = 0; j < other.get_array_size(); ++j) {
      result.at(i, j) = current - other.element(i, j);
    }
  }
  return result;
}
// _____________________________________________________________________________
DataField operator*(double current, DataField other) {  // friend
  DataField result(other.get_array_count(), other.get_array_size());
  for (size_t i = 0; i < other.get_array_count(); ++i) {
    for (size_t j = 0; j < other.get_array_size(); ++j) {
      result.at(i, j) = current * other.element(i, j);
    }
  }
  return result;
}
// _____________________________________________________________________________
DataField operator/(double current, DataField other) {  // friend
  DataField result(other.get_array_count(), other.get_array_size());
  for (size_t i = 0; i < other.get_array_count(); ++i) {
    for (size_t j = 0; j < other.get_array_size(); ++j) {
      result.at(i, j) = current / other.element(i, j);
    }
  }
  return result;
}
// _____________________________________________________________________________
// _____________________________________________________________________________
