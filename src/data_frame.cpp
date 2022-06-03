// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "data_frame.hpp"  // NOLINT
// _____________________________________________________________________________
DataFrame::DataFrame() {
  //
}
// _____________________________________________________________________________
DataFrame::DataFrame(const Properties& properties) {
  //
}
// _____________________________________________________________________________
DataFrame::DataFrame(const DataFrame& other) {
  //
}
// _____________________________________________________________________________
DataFrame::~DataFrame() {
  //
}
// _____________________________________________________________________________
DataFrame& DataFrame::operator=(const DataFrame& other) {
  std::cerr << "DataFrame::operator=():\"";
  std::cerr << "This is a base class function, which should never be called.\"";
  std::cerr << std::endl;
  return *this;
}
// _____________________________________________________________________________
DataFrame& DataFrame::operator+=(const DataFrame& other) {
  std::cerr << "DataFrame::operator+=():\"";
  std::cerr << "This is a base class function, which should never be called.\"";
  std::cerr << std::endl;
  return *this;
}
// _____________________________________________________________________________
DataFrame& DataFrame::operator-=(const DataFrame& other) {
  std::cerr << "DataFrame::operator-=():\"";
  std::cerr << "This is a base class function, which should never be called.\"";
  std::cerr << std::endl;
  return *this;
}
// _____________________________________________________________________________
DataFrame& DataFrame::operator*=(const DataFrame& other) {
  std::cerr << "DataFrame::operator*=():\"";
  std::cerr << "This is a base class function, which should never be called.\"";
  std::cerr << std::endl;
  return *this;
}
// _____________________________________________________________________________
DataFrame& DataFrame::operator/=(const DataFrame& other) {
  std::cerr << "DataFrame::operator/=():\"";
  std::cerr << "This is a base class function, which should never be called.\"";
  std::cerr << std::endl;
  return *this;
}
// _____________________________________________________________________________
DataFrame& DataFrame::operator*=(const double other) {
  std::cerr << "DataFrame::operator*=():\"";
  std::cerr << "This is a base class function, which should never be called.\"";
  std::cerr << std::endl;
  return *this;
}
// _____________________________________________________________________________
DataFrame DataFrame::operator+(const DataFrame& other) {
  std::cerr << "DataFrame::operator+():\"";
  std::cerr << "This is a base class function, which should never be called.\"";
  std::cerr << std::endl;
  return *this;
}
// _____________________________________________________________________________
DataFrame DataFrame::operator-(const DataFrame& other) {
  std::cerr << "DataFrame::operator-():\"";
  std::cerr << "This is a base class function, which should never be called.\"";
  std::cerr << std::endl;
  return *this;
}
// _____________________________________________________________________________
DataFrame DataFrame::operator*(const DataFrame& other) {
  std::cerr << "DataFrame::operator*():\"";
  std::cerr << "This is a base class function, which should never be called.\"";
  std::cerr << std::endl;
  return *this;
}
// _____________________________________________________________________________
DataFrame DataFrame::operator/(const DataFrame& other) {
  std::cerr << "DataFrame::operator/():\"";
  std::cerr << "This is a base class function, which should never be called.\"";
  std::cerr << std::endl;
  return *this;
}
// _____________________________________________________________________________
DataFrame DataFrame::operator*(const double other) {
  std::cerr << "DataFrame::operator*():\"";
  std::cerr << "This is a base class function, which should never be called.\"";
  std::cerr << std::endl;
  return *this;
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
