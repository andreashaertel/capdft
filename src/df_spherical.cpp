// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2022 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "df_spherical.hpp"  // NOLINT
#include <cmath>
#include <fftw3.h>
// Class type declarations
template class DFSpherical<double>;
template class DFSpherical<fftw_complex>;
// _____________________________________________________________________________
template <typename T>
size_t DFSpherical<T>::size() const {
  return array_size;
}
// _____________________________________________________________________________
template <typename T>
T& DFSpherical<T>::at(size_t i) {
  if (i >= array_size) {
    std::cerr << "DFSpherical::at(): \"Error: Index out of range.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return data[i];
}
// _____________________________________________________________________________
template <typename T>
T& DFSpherical<T>::element(size_t i) const {
  if (i >= array_size) {
    std::cerr << "DFSpherical::element():";
    std::cerr << "\"Error: Index out of range.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return data[i];
}
// _____________________________________________________________________________
template <typename T>
T* DFSpherical<T>::array() {
  return data;
}
// _____________________________________________________________________________
template <typename T>
void DFSpherical<T>::set_all_elements_to(T value) {
  for (size_t i = 0; i < array_size; ++i) {
    this->at(i) = value;
  }
}
template <>
void DFSpherical<fftw_complex>::set_all_elements_to(fftw_complex value) {
  for (size_t i = 0; i < array_size; ++i) {
    this->at(i)[0] = value[0];
    this->at(i)[1] = value[1];
  }
}
// _____________________________________________________________________________
template <typename T>
DFSpherical<T>::DFSpherical(size_t array_size) : array_size(array_size) {
  data = new T[array_size];
  set_all_elements_to(0.);
}
template <>
DFSpherical<fftw_complex>::DFSpherical(size_t array_size)
  : array_size(array_size) {
  data = new fftw_complex[array_size];
  fftw_complex zero{0.,0.};
  set_all_elements_to(zero);
}
// _____________________________________________________________________________
template <typename T>
DFSpherical<T>::DFSpherical(const Properties& properties) {
  properties.get_property("grid count", &array_size);
  data = new T[array_size];
  set_all_elements_to(0.);
}
template <>
DFSpherical<fftw_complex>::DFSpherical(const Properties& properties) {
  properties.get_property("grid count", &array_size);
  data = new fftw_complex[array_size];
  fftw_complex zero{0.,0.};
  set_all_elements_to(zero);
}
// _____________________________________________________________________________
template <typename T>
DFSpherical<T>::DFSpherical(const DFSpherical<T>& other) {
  this->array_size = other.size();
  data = new T[array_size];
  for (size_t i = 0; i < array_size; ++i) {
    this->at(i) = other.element(i);
  }
}
template <>
DFSpherical<fftw_complex>::DFSpherical(
    const DFSpherical<fftw_complex>& other) {
  this->array_size = other.size();
  data = new fftw_complex[array_size];
  for (size_t i = 0; i < array_size; ++i) {
    this->at(i)[0] = other.element(i)[0];  // Real
    this->at(i)[1] = other.element(i)[1];  // Imaginary
  }
}
//// _____________________________________________________________________________
template <typename T>
DFSpherical<T>::~DFSpherical() {
  delete [] data;
}
// _____________________________________________________________________________
template <typename T>
bool DFSpherical<T>::same_size(
    const DFSpherical<T>& other) const {
  return (this->size() == other.size());
}
// _____________________________________________________________________________
template <typename T>
void DFSpherical<T>::print(
    std::ostream& outstream, std::streamsize stream_size) const {
  // Save previous output precision
  std::streamsize old_stream_size = outstream.precision();
  // Header
  outstream.precision(stream_size);
  outstream << "# [index] [data]" << std::endl;
  // Data
  for (size_t i = 0; i < array_size; ++i) {
    outstream << i << " " << this->element(i) << std::endl;
  }
  // Restore output precision
  outstream.precision(old_stream_size);
}
template <>
void DFSpherical<fftw_complex>::print(
    std::ostream& outstream, std::streamsize stream_size) const {
  // Save previous output precision
  std::streamsize old_stream_size = outstream.precision();
  // Header
  outstream.precision(stream_size);
  outstream << "# [index] [data]" << std::endl;
  // Data
  for (size_t i = 0; i < array_size; ++i) {
    outstream << i << " ";
    outstream << this->element(i)[0];  // real part
    outstream << "+i" << this->element(i)[1];  // imaginary part
    outstream << std::endl;
  }
  // Restore output precision
  outstream.precision(old_stream_size);
}
// _____________________________________________________________________________
template <typename T>
void DFSpherical<T>::print(std::ostream& outstream) const {
  print(outstream, std::numeric_limits<double>::max_digits10);
}
// _____________________________________________________________________________
template <typename T>
void DFSpherical<T>::print() const {
  print(std::cout);
}
// _____________________________________________________________________________
template <typename T>
DFSpherical<T>& DFSpherical<T>::operator=(
    const DFSpherical<T>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) = other.element(i);
    }
  } else {
    std::cerr << "DFSpherical::operator=():";
    std::cerr << " \"ERROR: Array sizes do not match. Cannot copy.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DFSpherical<fftw_complex>& DFSpherical<fftw_complex>::operator=(
    const DFSpherical<fftw_complex>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i)[0] = other.element(i)[0];  // Real
      this->at(i)[1] = other.element(i)[1];  // Imaginary
    }
  } else {
    std::cerr << "DFSpherical::operator=():";
    std::cerr << " \"ERROR: Array sizes do not match. Cannot copy.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
// _____________________________________________________________________________
template <typename T>
DFSpherical<T>& DFSpherical<T>::operator+=(
    const DFSpherical<T>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) += other.element(i);
    }
  } else {
    std::cerr << "DFSpherical::operator+=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot add.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DFSpherical<fftw_complex>& DFSpherical<fftw_complex>::operator+=(
    const DFSpherical<fftw_complex>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i)[0] += other.element(i)[0];  // Real
      this->at(i)[1] += other.element(i)[1];  // Imaginary
    }
  } else {
    std::cerr << "DFSpherical::operator+=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot add.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
// _____________________________________________________________________________
template <typename T>
DFSpherical<T>& DFSpherical<T>::operator-=(
    const DFSpherical<T>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) -= other.element(i);
    }
  } else {
    std::cerr << "DFSpherical::operator-=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DFSpherical<fftw_complex>& DFSpherical<fftw_complex>::operator-=(
    const DFSpherical<fftw_complex>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i)[0] -= other.element(i)[0];  // Real
      this->at(i)[1] -= other.element(i)[1];  // Imaginary
    }
  } else {
    std::cerr << "DFSpherical::operator-=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
// _____________________________________________________________________________
template <typename T>
DFSpherical<T>& DFSpherical<T>::operator*=(
    const DFSpherical<T>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) *= other.element(i);
    }
  } else {
    std::cerr << "DFSpherical::operator*=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot multiply.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DFSpherical<fftw_complex>& DFSpherical<fftw_complex>::operator*=(
    const DFSpherical<fftw_complex>& other) {
  if (this->same_size(other)) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i)[0] = this->element(i)[0] * other.element(i)[0] -
          this->element(i)[1] * other.element(i)[1];  // Real
      this->at(i)[1] = this->element(i)[0] * other.element(i)[1] +
          this->element(i)[1] * other.element(i)[0];  // Imaginary
    }
  } else {
    std::cerr << "DFSpherical::operator*=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot multiply.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
// _____________________________________________________________________________
template <typename T>
DFSpherical<T>& DFSpherical<T>::operator/=(
    const DFSpherical<T>& other) {
  if (this->size() == other.size()) {
    for (size_t i = 0; i < this->array_size; ++i) {
      this->at(i) /= other.element(i);
    }
  } else {
    std::cerr << "DFSpherical::operator/=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot divide.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
template <>
DFSpherical<fftw_complex>& DFSpherical<fftw_complex>::operator/=(
    const DFSpherical<fftw_complex>& other) {
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
    std::cerr << "DFSpherical::operator/=():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot divide.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return *this;
}
// _____________________________________________________________________________
template <typename T>
DFSpherical<T>& DFSpherical<T>::operator*=(const T other) {
  for (size_t i = 0; i < this->array_size; ++i) {
    this->at(i) *= other;
  }
  return *this;
}
template <>
DFSpherical<fftw_complex>& DFSpherical<fftw_complex>::operator*=(
    const fftw_complex other) {
  for (size_t i = 0; i < this->array_size; ++i) {
    this->at(i)[0] = this->element(i)[0] * other[0] -
        this->element(i)[1] * other[1];  // Real
    this->at(i)[0] = this->element(i)[0] * other[1] +
        this->element(i)[1] * other[0];  // Imaginary
  }
  return *this;
}
// _____________________________________________________________________________
template <typename T>
DFSpherical<T> DFSpherical<T>::operator+(
    const DFSpherical<T>& other) {
  DFSpherical<T> result(*this);
  if (this->same_size(other)) {
    result += other;
  } else {
    std::cerr << "DFSpherical::operator+():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot add.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return result;
}
// _____________________________________________________________________________
template <typename T>
DFSpherical<T> DFSpherical<T>::operator-(
    const DFSpherical<T>& other) {
  DFSpherical<T> result(*this);
  if (this->same_size(other)) {
    result -= other;
  } else {
    std::cerr << "DFSpherical::operator-():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot subtract.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return result;
}
// _____________________________________________________________________________
template <typename T>
DFSpherical<T> DFSpherical<T>::operator*(
    const DFSpherical<T>& other) {
  DFSpherical<T> result(*this);
  if (this->same_size(other)) {
    result *= other;
  } else {
    std::cerr << "DFSpherical::operator*():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot multiply.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return result;
}
// _____________________________________________________________________________
template <typename T>
DFSpherical<T> DFSpherical<T>::operator/(
    const DFSpherical<T>& other) {
  DFSpherical<T> result(*this);
  if (this->same_size(other)) {
    result /= other;
  } else {
    std::cerr << "DFSpherical::operator/():";
    std::cerr << " \"ERROR: Dimensions do not match. Cannot divide.\"";
    std::cerr << std::endl;
    exit(1);
  }
  return result;
}
// _____________________________________________________________________________
template <typename T>
DFSpherical<T> DFSpherical<T>::operator*(const T other) {
  DFSpherical<T> result(*this);
  result *= other;
  return result;
}
// _____________________________________________________________________________
DFSpherical<double> operator*(const double current,
    const DFSpherical<double>& other) {
  DFSpherical<double> result(other);
  result *= current;
  return result;
}
// _____________________________________________________________________________
DFSpherical<double> exp(const DFSpherical<double>& other) {
  DFSpherical<double> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i) = exp(other.element(i));
  }
  return result;
}
DFSpherical<fftw_complex> exp(const DFSpherical<fftw_complex>& other) {
  DFSpherical<fftw_complex> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i)[0] = exp(other.element(i)[0]) * cos(other.element(i)[1]);
    result.at(i)[1] = exp(other.element(i)[0]) * sin(other.element(i)[1]);
  }
  return result;
}
// _____________________________________________________________________________
DFSpherical<double> log_natural(const DFSpherical<double>& other) {
  DFSpherical<double> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i) = log(other.element(i));
  }
  return result;
}
DFSpherical<fftw_complex> log_natural(
    const DFSpherical<fftw_complex>& other) {
  double re{0.}, im{0.}, radius{0.}, angle{0.};
  DFSpherical<fftw_complex> result(other);
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
  std::cerr << "DFSpherical::log_natural():";
  std::cerr << " \"Warning: The complex logarithm is used.";
  std::cerr << "Keep in mind, that there are infinite solutions.\"";
  std::cerr << std::endl;
  return result;
}
// _____________________________________________________________________________
DFSpherical<double> abs(const DFSpherical<double>& other) {
  DFSpherical<double> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i) = fabs(result.at(i));
  }
  return result;
}
DFSpherical<fftw_complex> abs(const DFSpherical<fftw_complex>& other) {
  DFSpherical<fftw_complex> result(other);
  for (size_t i = 0; i < other.array_size; ++i) {
    result.at(i)[0] = sqrt(pow(result.at(i)[0], 2) + pow(result.at(i)[1], 2));
    result.at(i)[1] = 0.;
  }
  return result;
}
// _____________________________________________________________________________
double max(const DFSpherical<double>& other) {
  return *std::max_element(other.data, other.data + other.array_size);
}
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
