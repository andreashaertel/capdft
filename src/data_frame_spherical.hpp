// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2022 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_DATA_FRAME_SPHERICAL_HPP_
#define SRC_DATA_FRAME_SPHERICAL_HPP_
/** \file data_frame.hpp
 *  \brief This file contains the declarations of the DataFrameSpherical class,
 *  which is a child of DataFrame.
 *
 */
#include "src/properties.hpp"
#include "src/data_frame.hpp"
/** \brief Container class for general data (e.g. density profiles or functional
 *  derivatives) derived from DataFrame
 *
 *  The DataFrameSpherical is a child of the general data object of the capDFT
 *  project DataFrame, which can hold data profiles as density profiles,
 *  functional derivatives, or external potential fields.
 *  The DataFrameSpherical provides the specific implementation of the
 *  DataFrame class for the spherical geometry.
 *  It contains all virtual functions of the parent and 
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
  DataFrameSpherical(const DataFrameSpherical<T>& other);
  /** \brief Virtual destructor
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
   *  You can also modify an element in that way.
   *
   *  \return The array element reference at position i.
   *
   */
  T& at(size_t i);
  /** \brief Access the array element at position i like at(),
   *  but without modifying it.
   *
   *  \return The array element reference at position i.
   *
   */
  T& element(size_t i) const;
  /** \brief Test for DataFrameSpherical other having the same size as this. 
   *
   *  \param other The DataFrameSpherical thats size is compared to this DataFrameSpherical. 
   *
   *  \return True, if the other DataFrameSpherical has the same size as this one, 
   *          False otherwise. 
   *
   */
  virtual bool same_size(const DataFrameSpherical<T>& other) const;
  /** \brief The = operator copies the content of other into this. */
  virtual DataFrameSpherical<T>& operator=(const DataFrameSpherical<T>& other);
  /** \brief The += operator adds the content of other to this and then returns
   *         this. 
   */
  virtual DataFrameSpherical<T>& operator+=(const DataFrameSpherical<T>& other);
  /** \brief The -= operator subtracts the content of other from this and then
   *         returns this. 
   */
  virtual DataFrameSpherical<T>& operator-=(const DataFrameSpherical<T>& other);
  /** \brief The *= operator multiplies the content of other to this and then
   *         returns this. 
   */
  virtual DataFrameSpherical<T>& operator*=(const DataFrameSpherical<T>& other);
  /** \brief The /= operator divides this by the content of other and then 
   *         returns this. 
   */
  virtual DataFrameSpherical<T>& operator/=(const DataFrameSpherical<T>& other);
  /** \brief The *= operator for scalar multiplication multiplies all entries of
   *         this with the scalar value other and then returns this. 
   */
  virtual DataFrameSpherical<T>& operator*=(const double other);
  /** \brief The + operator adds the content of this and other and returns the 
   *         result. 
   *
   *  Use the following code: 
   *
   *      DataFrameSpherical& operator+(const DataFrameSpherical& other) {
   *        // Check for correct sizes
   *        if (! this->same_size(other))
   *        throw &bad_size_error;
   *        // Use existing copy Constructor and += operators.
   *        DataFrameSpherical result(*this);
   *        result += other;
   *        return result;
   *      };
   *
   */
  virtual DataFrameSpherical<T>& operator+(const DataFrameSpherical<T>& other)
      = 0;
  /** \brief The - operator subtracts the content of this and other and returns
   *         the result. 
   *
   *  Use the following code: 
   *
   *      DataFrameSpherical& DataFrameSpherical::operator-(const DataFrameSpherical& other) {
   *        // Check for correct sizes
   *        if (! this->same_size(other))
   *        throw &bad_size_error;
   *        // Use existing copy Constructor and -= operators.
   *        DataFrameSpherical result(*this);
   *        result -= other;
   *        return result;
   *      }
   *
   */
  virtual DataFrameSpherical<T>& operator-(const DataFrameSpherical<T>& other)
      = 0;
  /** \brief The * operator multiplies the content of this and other and returns
   *         the result. 
   *
   *  Use the following code: 
   *
   *      DataFrameSpherical& DataFrameSpherical::operator*(const DataFrameSpherical& other) {
   *        // Check for correct sizes
   *        if (! this->same_size(other))
   *        throw &bad_size_error;
   *        // Use existing copy Constructor and *= operators.
   *        DataFrameSpherical result(*this);
   *        result *= other;
   *        return result;
   *      }
   *
   */
  virtual DataFrameSpherical<T>& operator*(const DataFrameSpherical<T>& other)
      = 0;
  /** \brief The / operator divides the content of this and other and returns 
   *         the result. 
   *
   *  Use the following code: 
   *
   *      DataFrameSpherical& DataFrameSpherical::operator/(const DataFrameSpherical& other) {
   *      // Check for correct sizes
   *      if (! this->same_size(other))
   *      throw &bad_size_error;
   *      // Use existing copy Constructor and /= operators.
   *      DataFrameSpherical result(*this);
   *      result /= other;
   *      return result;
   *    }
   *
   */
  virtual DataFrameSpherical<T>& operator/(const DataFrameSpherical<T>& other)
      = 0;
  /** \brief The * operator multiplies the content of this and the scalar other
   *         and returns the result. 
   *
   *  Use the following code: 
   *
   *      DataFrameSpherical& DataFrameSpherical::operator*(const double other) {
   *        // Use existing copy Constructor and *= for double operators.
   *        DataFrameSpherical result(*this);
   *        result *= other;
   *        return result;
   *      }
   *
   */
  virtual DataFrameSpherical<T>& operator*(const double other) = 0;
  /** \brief Exponential function applied to all elements of DataFrameSpherical. 
   *
   */
  virtual DataFrameSpherical<T>& exp() = 0;
  /** \brief Natural logarithm function applied to all elements of
   *  DataFrameSpherical. 
   *
   */
  virtual DataFrameSpherical<T>& log_natural() = 0;
  /** \brief The * operator for (double) * DataFrameSpherical. 
   *
   *  Use the following code: 
   *
   *      friend DataFrameSpherical& operator*(const double current,
   *          const DataFrameSpherical& other) {
   *        // Use existing copy Constructor and *= for double operators.
   *        DataFrameSpherical result(other);
   *        result *= current;
   *        return result;
   *      }
   *
   */
  template <typename U>
  friend DataFrameSpherical<U>& operator*(
      const double current, const DataFrameSpherical<U>& other);
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
