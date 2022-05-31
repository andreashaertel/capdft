// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2022 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_DATA_FRAME_SPHERICAL_HPP_
#define SRC_DATA_FRAME_SPHERICAL_HPP_
/** \file data_frame.hpp
 *  \brief This file contains the declarations of the DataFrameSpherical class,
 *  which is a child class of DataFrame.
 *
 */
#include "src/properties.hpp"
#include "src/data_frame.hpp"
/** \brief Container class for general data (e.g. density profiles or functional
 *         derivatives) derived from DataFrame
 *
 *  The DataFrameSpherical is a child of the general data object of the capDFT
 *  project DataFrame, which can hold data profiles as density profiles,
 *  functional derivatives, or external potential fields.
 *  The DataFrameSpherical provides the specific implementation of the
 *  DataFrame class for the spherical geometry.
 *  It contains all virtual functions of the parent and some additional quality
 *  of life functions. The use of templates makes the extension to other vaiable
 *  types like "int" easier.
 *
 */
template <typename T>
class DataFrameSpherical : public DataFrame {
 public:
  /** \brief Constructors
   *
   *  A DataFrameSpherical has two Constructors: 
   *  One to construct a new implementation, and a second one to create a copy 
   *  of an existing one. To create a new implementation, parameters are handed
   *  over to the constructor as Properties.
   *
   *  The copy Constructor, however, just creates a copy of an existing
   *  DataFrameSpherical.
   *
   */
  explicit DataFrameSpherical(const Properties& properties);
  DataFrameSpherical(const DataFrame& other);
  /** \brief Destructor
   *
   */
  virtual ~DataFrameSpherical();
  /** \brief Return internal array size. 
   *
   *  \return The size of the internal array as size_t. 
   *
   */
  size_t size() const;
  /** \brief Access the array element at position i.
   *         You can also modify an element in that way.
   *
   *  Writing: my_object.at(i) = 5.;
   *  Reading: my_var = my_object.at(i);
   *
   *  \return The array element reference at position i.
   *
   */
  T& at(size_t i);
  /** \brief Access the array element at position i like at(),
   *         but without modifying it.
   *
   *  \return The array element reference at position i.
   *
   */
  T& element(size_t i) const;
  /** \brief Test if the DataFrameSpherical object other has the same size as
   *         this object. 
   *
   *  \param other The DataFrameSpherical thats size is compared to this
   *         DataFrameSpherical. 
   *
   *  \return true, if the other DataFrameSpherical has the same size as this
   *          one, false otherwise. 
   *
   */
  virtual bool same_size(const DataFrameSpherical<T>& other) const;
  /** \brief The = operator copies the content of other into this object.
   *
   *  \return DataFrameSpherical reference of this object.
   *
   **/
  virtual DataFrameSpherical<T>& operator=(const DataFrame& other);
  /** \brief The += operator adds the content of other to this and then returns
   *         this object reference. 
   *
   *  \return DataFrameSpherical reference of this object.
   *
   */
  virtual DataFrameSpherical<T>& operator+=(const DataFrame& other);
  /** \brief The -= operator subtracts the content of other from this and then
   *         returns this. 
   *
   *  \return DataFrameSpherical reference of this object.
   *
   */
  virtual DataFrameSpherical<T>& operator-=(const DataFrame& other);
  /** \brief The *= operator multiplies the content of other to this and then
   *         returns this. 
   *
   *  \return DataFrameSpherical reference of this object.
   *
   */
  virtual DataFrameSpherical<T>& operator*=(const DataFrame& other);
  /** \brief The /= operator divides this by the content of other and then 
   *         returns this. 
   *
   *  \return DataFrameSpherical reference of this object.
   *
   */
  virtual DataFrameSpherical<T>& operator/=(const DataFrame& other);
  /** \brief The *= operator for scalar multiplication multiplies all entries of
   *         this with the scalar value other and then returns this. 
   *
   *  \return DataFrameSpherical reference of this object.
   *
   */
  virtual DataFrameSpherical<T>& operator*=(const double other);
  /** \brief The + operator adds the content of this object and the
   *         DataFrame other and returns the result. 
   *
   *  \return DataFrame object which is the sum of this object and
   *          and the DataFrameSpherical other.
   *
   */
  virtual DataFrame operator+(const DataFrame& other);
  /** \brief The - operator subtracts the content of this object and the
   *         DataFrame other and returns the result. 
   *
   *  \return DataFrame object which is the difference of this object
   *          and the DataFrame other.
   *
   */
  virtual DataFrame operator-(const DataFrame& other);
  /** \brief The * operator multiplies the content of this object and the
   *         DataFrame other and returns the result. 
   *
   *  \return DataFrame object which is the product of this object and
   *          the DataFrame other.
   *
   */
  virtual DataFrame operator*(const DataFrame& other);
  /** \brief The / operator divides the content of this object and the DataFrame
   *         other and returns the result. 
   *
   *  \return DataFrame object which is the quotient of this object and the
   *          DataFrame other.
   *
   */
  virtual DataFrame operator/(const DataFrame& other);
  /** \brief The * operator multiplies the content of this object and the scalar
   *         other and returns the result. 
   *
   *  \return DataFrame object which is the product this object and
   *          the double other.
   *
   */
  virtual DataFrame operator*(const T other);
  /** \brief The * operator for (<typename T>) * DataFrameSpherical. 
   *
   *  \return DataFrameSpherical object which is the product of a the variable
   *          current and the DataFrameSpherical other.
   *
   */
  template <typename U>
  friend DataFrameSpherical<U> operator*(const U current,
      const DataFrameSpherical<U>& other);
  /** \brief Exponential function that returns a DataFrameSpherical where all
   *         elements of DataFrameSpherical other were exponentiated.
   *
   *  \return DataFrameSpherical with exponentiated values of the
   *          DataFrameSpherical other.
   *
   */
  template <typename V>
  friend DataFrameSpherical<V> exp(DataFrameSpherical<V>& other);
  /** \brief Natural logarithm that returns a DataFrameSpherical where all
   *         elements of DataFrameSpherical other were logarithmized.
   *
   *  \return DataFrameSpherical with logarithmized values of the
   *          DataFrameSpherical other.
   *
   */
  template <typename W>
  friend DataFrameSpherical<W> log_natural(DataFrameSpherical<W>& other);
  /** \brief std::exception BadSizeException.
   */
  class BadSizeException : public std::exception {
   public:
    /** \brief Overwrite the exception information function what().
     *
     *  \return the text "Size of DataFrameSphericals does not match.".
     */
    virtual const char* what(void) const throw() {
      return "Size of DataFrameSphericals does not match.";
    }
  }
  /** \brief Exception BadSizeException bad_size_error.
   *
   *  The exception is thrown if two DataFrameSphericals do not match in size. 
   */
  bad_size_error;

 protected:
  //

 private:
 /** \brief Size of the internal array
  */
  size_t array_size;
 /** \brief Array holding the data of size array_size
  */
  T* data;
};
#endif  // SRC_DATA_FRAME_SPHERICAL_HPP_
