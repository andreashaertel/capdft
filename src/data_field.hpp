// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_DATA_FIELD_HPP_
#define SRC_DATA_FIELD_HPP_
/** \file data_field.cpp
 *  \brief This file contains the declarations of the DataField class
 *
 */
#include <cstddef>
#include <vector>
#include <string>
/** \brief Container class for arrays (e.g. densities, functional derivatives)
 *
 */
class DataField {
 public:
  /** \brief Constructors
   *
   * There are several types of constructors:
   * 1. A copy constructor
   * 2. An empty constructor
   * 3. One that does nnot set a bin width (standard)
   * 4. One that also sets a bin width
   *
   */
  DataField();
  DataField(size_t array_count, size_t array_size);
  DataField(size_t array_count, size_t array_size, double bin_width);
  DataField(const DataField& other);
  /** \brief Destructor
   *
   */
  ~DataField();
  /** \brief Access to one array element j of a certain array i (read/write)
   *
   */
  double& at(size_t i, size_t j);
  /** \brief Returns one array element j of a certain array i (read only)
   *
   */
  double& element(size_t i, size_t j) const;
  /** \brief Access to one array i
   *
   */
  double* array(size_t i);
  /** \brief Returns outer iterators of the vector
   *
   */
  std::vector<double*>::iterator begin();
  std::vector<double*>::iterator end();
  /** \brief Free all allocated memory
   *
   */
  void clear();
  /** \brief Allocate required memory. Used by most constructors.
   *
   */
  void create(size_t array_count, size_t array_size);
  /** \brief Do clear() and create()
   *
   */
  void recreate(size_t array_count, size_t array_size);
  /** \brief Returns the number of arrays
   *
   */
  size_t get_array_count() const;
  /** \brief Returns the size of the arrays
   *
   */
  size_t get_array_size() const;
  /** \brief Returns bin width
   *
   */
  double get_bin_width() const;
  /** \brief Sets bin width
   *
   */
  void set_bin_width(double bw);
  /** \brief Print arrays to specified output with custom or standard precision
   *
   */
  void print(std::ostream& outstream, std::streamsize stream_size) const;
  void print(std::ostream& outstream) const;
  /** \brief Print all values to standard output
   *
   */
  void print() const;
  /** \brief Write a certain number to all arrays
   *
   */
  void set_all_elements_to(const double number);
  /** \brief Write zeros to all arrays
   *
   */
  void zeros();
  /** \brief The = operator copies the whole object
   *
   *  List of allowed operator actions:
   *  1.) DataField = DataField
   *  2.) DataField += DataField
   *  3.) DataField -= DataField
   *  4.) DataField *= <number>
   *  5.) DataField /= <number>
   *  6.) DataField + DataField
   *  7.) DataField - DataField
   *  8.) DataField + <number>
   *  9.) DataField - <number>
   *  10.) DataField * <number>
   *  11.) DataField / <number>
   *
   */
  DataField& operator=(DataField other);
  /** \brief The += operator
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField& operator+=(DataField other);
  /** \brief The -= operator
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField& operator-=(DataField other);
  /** \brief The *= operator
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField& operator*=(double other);
  /** \brief The /= operator
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField& operator/=(double other);
  /** \brief The + operator adds two DataFields elementwise
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField operator+(DataField other);
  /** \brief The - operator subtracts two DataFields elementwise
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField operator-(DataField other);
  /** \brief The + operator adds a constant value to every entry
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField operator+(double other);  // for commuted case see bottom
  /** \brief The - operator subtracts a constant value from every entry
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField operator-(double other);  // for commuted case see bottom
  /** \brief The * operator multiplies every entry by a constant factor
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField operator*(double other);  // for commuted case see bottom
  /** \brief The / operator divides every entry by a constant value
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField operator/(double other);  // for commuted case see bottom
  /** \brief Addition operator has to be defined both ways for commutativity
   *
   *  See list of all operator actions at operator=().
   *
   */
  friend DataField operator+(double current, DataField other);
  /** \brief Subtraction operator has to be defined both ways for commutativity
   *
   *  See list of all operator actions at operator=().
   *
   */
  friend DataField operator-(double current, DataField other);
  /** \brief Multiplication operator has to be defined
   *
   *  See list of all operator actions at operator=().
   *
   */
  friend DataField operator*(double current, DataField other);
  /** \brief Division operator has to be defined both ways
   *
   *  See list of all operator actions at operator=().
   *
   */
  friend DataField operator/(double current, DataField other);

 private:
  size_t array_count;
  size_t array_size;
  double bin_width;
  std::vector<double*> arrays;

 protected:
};
#endif  // SRC_DATA_FIELD_HPP_
