// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2022 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_DATA_FRAME_HPP_
#define SRC_DATA_FRAME_HPP_
/** \file data_frame.hpp
 *  \brief This file contains the declarations of the DataFrame class,
 *  which is a certain data container class.
 *
 */
#include <algorithm>
#include <string>
#include <vector>
/** \brief Container class for general data (e.g. density profiles or functional
 *         derivatives)
 *
 *  The DataFrame is the data class of capDFT, which can hold data like density
 *  profiles, functional derivatives, or external potential fields.
 *  There are two template arguments. The first one determines the
 *  dimensionality of your DataFrame. The second one determines the underlying
 *  type of the data to be stored (i.e. "double", "fftw_complex").
 *
 */
template <size_t dim, typename T>
class DataFrame {
 public:
  /** \brief Constructors
   *
   *  A DataFrame has several constructors:
   *  The first one is the standard constructor, which creates an empty
   *  DataFrame.
   *  The second one creates a 1D-DataFrame with array_size bins.
   *  The third one creates a data frame of arbitrary dimensionality. The bin
   *  in each dimension is given by the supplied vector
   *  (e.g. the input {3,3,3} creates a 3 by 3 by 3 grid).
   *  The fourth one is the copy constructor, which creates a copy of an
   *  existing DataFrame.
   *
   */
  DataFrame();
  explicit DataFrame(size_t array_size);
  explicit DataFrame(std::vector<size_t> array_dimensions);
  DataFrame(const DataFrame<dim, T>& other);
  /** \brief Destructor
   *
   */
  ~DataFrame();
  /** \brief Return internal array size. 
   *
   *  \return The size of the internal array as size_t. 
   *
   */
  size_t size() const;
  /** \brief Return the grid count of each dimension.
   *
   *  \return The grid count of each dimension as vector.
   *
   */
  std::vector<size_t> size_dim() const;
  /** \brief Access array elements.
   *  
   *  at(i) works in for all dimensions, since the data is always a 1D array.
   *  at(i,j) and at(i,j,k) only work in 2 and 3 dimensions respectively.
   *  You can read and write with at().
   *
   *  Writing: my_object.at(i) = 5.;
   *  Reading: my_var = my_object.at(i);
   *
   *  \return The array element reference at position i (or (i,j), or (i,j,k)).
   *
   */
  T& at(size_t i);
  T& at(size_t i, size_t j);
  T& at(size_t i, size_t j, size_t k);
  /** \brief Access array elements like at(), but without the ability to modify
   *         them.
   *
   *  \return The array element reference at position i (or (i,j), or (i,j,k)).
   *
   */
  T& element(size_t i) const;
  T& element(size_t i, size_t j) const;
  T& element(size_t i, size_t j, size_t k) const;
  /** \brief Return element i with precision stream_size as string.
   *
   *  This function is used in the print()-function.
   *
   *  \return The array element as string.
   *
   */
  std::string element_string(size_t i, std::streamsize stream_size) const;
  /** \brief Function that returns the reference to the data array pointer.
   *
   */
  T*& array();
  /** \brief Function that sets all elements in the DataFrame object to
   *         a certain value.
   *
   */
  void set_all_elements_to(T value);
  /** \brief Function that sets all elements in the DataFrame object to 0.
   *
   */
  void zero();
  /** \brief Test if the DataFrame object other has the same size as
   *         this object. 
   *
   *  \param other The DataFrame thats size is compared to this
   *         DataFrame. 
   *
   *  \return true, if the other DataFrame has the same size as this
   *          one, false otherwise. 
   *
   */
  bool same_size(const DataFrame<dim, T>& other) const;
  /** \brief Print arrays to specified output with custom or standard precision
   *
   */
  void print(std::ostream& outstream, std::streamsize stream_size) const;
  void print(std::ostream& outstream) const;
  /** \brief Print all values to standard output
   *
   */
  void print() const;
  /** \brief The = operator copies the content of other into this object.
   *
   *  \return DataFrame reference of this object.
   *
   **/
  DataFrame<dim, T>& operator=(const DataFrame<dim, T>& other);
  /** \brief The += operator adds the content of other to this and then returns
   *         this object reference. 
   *
   *  \return DataFrame reference of this object.
   *
   */
  DataFrame<dim, T>& operator+=(const DataFrame<dim, T>& other);
  /** \brief The -= operator subtracts the content of other from this and then
   *         returns this. 
   *
   *  \return DataFrame reference of this object.
   *
   */
  DataFrame<dim, T>& operator-=(const DataFrame<dim, T>& other);
  /** \brief The *= operator multiplies the content of other to this and then
   *         returns this. 
   *
   *  \return DataFrame reference of this object.
   *
   */
  DataFrame<dim, T>& operator*=(const DataFrame<dim, T>& other);
  /** \brief The /= operator divides this by the content of other and then 
   *         returns this. 
   *
   *  \return DataFrame reference of this object.
   *
   */
  DataFrame<dim, T>& operator/=(const DataFrame<dim, T>& other);
  /** \brief The *= operator for scalar multiplication multiplies all entries of
   *         this with the scalar value other and then returns this. 
   *
   *  \return DataFrame reference of this object.
   *
   */
  DataFrame<dim, T>& operator*=(const double other);
  /** \brief The + operator adds the content of this object and the
   *         DataFrame other and returns the result. 
   *
   *  \return DataFrame object which is the sum of this object and
   *          and the DataFrame other.
   *
   */
  DataFrame<dim, T> operator+(const DataFrame<dim, T>& other);
  /** \brief The - operator subtracts the content of this object and the
   *         DataFrame other and returns the result. 
   *
   *  \return DataFrame object which is the difference of this object
   *          and the DataFrame other.
   *
   */
  DataFrame<dim, T> operator-(const DataFrame<dim, T>& other);
  /** \brief The * operator multiplies the content of this object and the
   *         DataFrameother and returns the result. 
   *
   *  \return DataFrameobject which is the product of this object and
   *          the DataFrameother.
   *
   */
  DataFrame<dim, T> operator*(const DataFrame<dim, T>& other);
  /** \brief The / operator divides the content of this object and the DataFrame
   *         other and returns the result. 
   *
   *  \return DataFrameobject which is the quotient of this object and the
   *          DataFrameother.
   *
   */
  DataFrame<dim, T> operator/(const DataFrame<dim, T>& other);
  /** \brief The * operator multiplies the content of this object and the scalar
   *         other and returns the result. 
   *
   *  \return DataFrameobject which is the product this object and
   *          the the template value T other.
   *
   */
  DataFrame<dim, T> operator*(const double other);
  /** \brief The * operator for (<typename T>) * DataFrame. 
   *
   *  \return DataFrame object which is the product of a the variable
   *          current and the DataFrame other.
   *
   */
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wnon-template-friend"
  friend DataFrame<dim, T> operator*(const double current,
      const DataFrame<dim, T>& other);
  /** \brief Exponential function that returns a DataFrame where all
   *         elements of DataFrame other were exponentiated.
   *
   *  \return DataFrame with exponentiated values of the
   *          DataFrame other.
   *
   */
  friend DataFrame<dim, T> exp(const DataFrame<dim, T>& other);
  /** \brief Natural logarithm that returns a DataFrame where all
   *         elements of DataFrame other were logarithmized.
   *
   *  \return DataFrame with logarithmized values of the
   *          DataFrame other.
   *
   */
  friend DataFrame<dim, T> log_natural(const DataFrame<dim, T>& other);
  /** \brief Apply the modulus (absolute value) to all elements.
   *
   */
  friend DataFrame<dim, T> abs(const DataFrame<dim, T>& other);
  /** \brief Return the largest element of the array.
   *         Makes only sense for data of type double.
   *
   */
  friend double max(const DataFrame<dim, double>& other);
  #pragma GCC diagnostic pop
  /** \brief std::exception BadSizeException.
   */
  class BadSizeException : public std::exception {
   public:
    /** \brief Overwrite the exception information function what().
     *
     *  \return the text "Size of DataFrames does not match.".
     */
    virtual const char* what(void) const throw() {
      return "Size of DataFrames does not match.";
    }
  }
  /** \brief Exception BadSizeException bad_size_error.
   *
   *  The exception is thrown if two DataFrames do not match in size. 
   */
  bad_size_error;

 protected:
  //

 private:
 /** \brief Total size of the internal array
  *
  */
  size_t array_size;
 /** \brief Vector that holds the size in every dimension separately
  *
  */
  std::vector<size_t> array_dimensions;
 /** \brief Array holding the data of size array_size
  */
  T* data;
  /** \brief Function that converts the total index into dim coordinates
   *
   */
  std::vector<size_t> index_to_coordinates(size_t index) const;
  /** \brief Function that converts dim coordinates to the total index.
   *
   */
  size_t coordinates_to_index(std::vector<size_t> coordinates) const;
};
#endif  // SRC_DATA_FRAME_HPP_
