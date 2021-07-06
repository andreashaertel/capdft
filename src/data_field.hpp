// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_DATA_FIELD_HPP_
#define SRC_DATA_FIELD_HPP_
/** \file data_field.hpp
 *  \brief This file contains the declarations of the DataField class
 *
 */
#include <cstddef>
#include <fftw3.h>
#include <vector>
#include <string>
/** \brief Container class for arrays (e.g. densities, functional derivatives)
 *
 */
template <typename T>
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
  DataField(const DataField<T>& other);
  /** \brief Destructor
   *
   */
  ~DataField();
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
  /** \brief Access to one array element j of a certain array i (read/write)
   *
   */
  T& at(size_t i, size_t j);
  /** \brief Returns one array element j of a certain array i (read only)
   *
   */
  T& element(size_t i, size_t j) const;
  /** \brief Pointer access to one array i
   *
   */
  T* array(size_t i);
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
   *  4.) DataField *= DataField
   *  5.) DataField /= DataField
   *  6.) DataField += <number>
   *  7.) DataField -= <number>
   *  8.) DataField *= <number>
   *  9.) DataField /= <number>
   *  10.) DataField + DataField
   *  11.) DataField - DataField
   *  12.) DataField * DataField
   *  13.) DataField / DataField
   *  14.) DataField + <number>
   *  15.) DataField - <number>
   *  16.) DataField * <number>
   *  17.) DataField / <number>
   *
   */
  DataField<T>& operator=(DataField<T> other);
  /** \brief The += operator
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField<T>& operator+=(DataField<T> other);
  /** \brief The -= operator
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField<T>& operator-=(DataField<T> other);
  /** \brief The *= operator
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField<T>& operator*=(DataField<T> other);
  /** \brief The /= operator
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField<T>& operator/=(DataField<T> other);
  /** \brief The += operator
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField<T>& operator+=(double other);
  /** \brief The -= operator
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField<T>& operator-=(double other);
  /** \brief The *= operator
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField<T>& operator*=(double other);
  /** \brief The /= operator
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField<T>& operator/=(double other);
  /** \brief The + operator adds two DataFields elementwise
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField<T> operator+(DataField<T> other);
  /** \brief The - operator subtracts two DataFields elementwise
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField<T> operator-(DataField<T> other);
  /** \brief The * operator multiplies two DataFields elementwise
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField<T> operator*(DataField<T> other);
  /** \brief The / operator divides a DataField by another DataField elementwise
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField<T> operator/(DataField<T> other);
  /** \brief The + operator adds a constant value to every entry
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField<T> operator+(double other);  // for commuted case see bottom
  /** \brief The - operator subtracts a constant value from every entry
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField<T> operator-(double other);  // for commuted case see bottom
  /** \brief The * operator multiplies every entry by a constant factor
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField<T> operator*(double other);  // for commuted case see bottom
  /** \brief The / operator divides every entry by a constant value
   *
   *  See list of all operator actions at operator=().
   *
   */
  DataField<T> operator/(double other);  // for commuted case see bottom
  /** \brief Addition operator has to be defined both ways for commutativity
   *
   *  See list of all operator actions at operator=().
   *
   */
  template <typename U>
  friend DataField<U> operator+(double current, DataField<U> other);
  /** \brief Subtraction operator has to be defined both ways for commutativity
   *
   *  See list of all operator actions at operator=().
   *
   */
  template <typename V>
  friend DataField<V> operator-(double current, DataField<V> other);
  /** \brief Multiplication operator has to be defined
   *
   *  See list of all operator actions at operator=().
   *
   */
  template <typename W>
  friend DataField<W> operator*(double current, DataField<W> other);
  /** \brief Division operator has to be defined both ways
   *
   *  See list of all operator actions at operator=().
   *
   */
  template <typename X>
  friend DataField<X> operator/(double current, DataField<X> other);

 protected:
  size_t array_count;
  size_t array_size;
  double bin_width;
  bool allocated_memory;
  std::vector<T*> arrays;
};
#endif  // SRC_DATA_FIELD_HPP_
