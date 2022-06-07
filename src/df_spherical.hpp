// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2022 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_DF_SPHERICAL_HPP_
#define SRC_DF_SPHERICAL_HPP_
/** \file data_frame.hpp
 *  \brief This file contains the declarations of the DFSpherical class,
 *  which is a certain data container class.
 *
 */
#include "properties.hpp"  // NOLINT
/** \brief Container class for general data (e.g. density profiles or functional
 *         derivatives)
 *
 *  The DFSpherical is a data class of capDFT for spherical geometries,
 *  which can hold data like density profiles,
 *  functional derivatives, or external potential fields.
 *  The DFSpherical provides the specific implementation of the
 *  DFSpherical class for the spherical geometry.
 *  It contains all virtual functions of the parent and some additional quality
 *  of life functions. The use of templates makes the extension to other vaiable
 *  types like "int" easier.
 *
 */
template <typename T>
class DFSpherical {
 public:
  /** \brief Constructors
   *
   *  A DFSpherical has two Constructors: 
   *  One to construct a new implementation, and a second one to create a copy 
   *  of an existing one. To create a new implementation, parameters are handed
   *  over to the constructor as Properties.
   *
   *  The copy Constructor, however, just creates a copy of an existing
   *  DFSpherical.
   *
   */
  explicit DFSpherical(const Properties& properties);
  DFSpherical(const DFSpherical<T>& other);
  /** \brief Destructor
   *
   */
  ~DFSpherical();
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
  /** \brief Test if the DFSpherical object other has the same size as
   *         this object. 
   *
   *  \param other The DFSpherical thats size is compared to this
   *         DFSpherical. 
   *
   *  \return true, if the other DFSpherical has the same size as this
   *          one, false otherwise. 
   *
   */
  bool same_size(const DFSpherical<T>& other) const;
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
   *  \return DFSpherical reference of this object.
   *
   **/
  DFSpherical<T>& operator=(const DFSpherical<T>& other);
  /** \brief The += operator adds the content of other to this and then returns
   *         this object reference. 
   *
   *  \return DFSpherical reference of this object.
   *
   */
  DFSpherical<T>& operator+=(const DFSpherical<T>& other);
  /** \brief The -= operator subtracts the content of other from this and then
   *         returns this. 
   *
   *  \return DFSpherical reference of this object.
   *
   */
  DFSpherical<T>& operator-=(const DFSpherical<T>& other);
  /** \brief The *= operator multiplies the content of other to this and then
   *         returns this. 
   *
   *  \return DFSpherical reference of this object.
   *
   */
  DFSpherical<T>& operator*=(const DFSpherical<T>& other);
  /** \brief The /= operator divides this by the content of other and then 
   *         returns this. 
   *
   *  \return DFSpherical reference of this object.
   *
   */
  DFSpherical<T>& operator/=(const DFSpherical<T>& other);
  /** \brief The *= operator for scalar multiplication multiplies all entries of
   *         this with the scalar value other and then returns this. 
   *
   *  \return DFSpherical reference of this object.
   *
   */
  DFSpherical<T>& operator*=(const T other);
  /** \brief The + operator adds the content of this object and the
   *         DFSpherical other and returns the result. 
   *
   *  \return DFSpherical object which is the sum of this object and
   *          and the DFSpherical other.
   *
   */
  DFSpherical<T> operator+(const DFSpherical<T>& other);
  /** \brief The - operator subtracts the content of this object and the
   *         DFSpherical other and returns the result. 
   *
   *  \return DFSpherical object which is the difference of this object
   *          and the DFSpherical other.
   *
   */
  DFSpherical<T> operator-(const DFSpherical<T>& other);
  /** \brief The * operator multiplies the content of this object and the
   *         DFSphericalother and returns the result. 
   *
   *  \return DFSphericalobject which is the product of this object and
   *          the DFSphericalother.
   *
   */
  DFSpherical<T> operator*(const DFSpherical<T>& other);
  /** \brief The / operator divides the content of this object and the DataFrame
   *         other and returns the result. 
   *
   *  \return DFSphericalobject which is the quotient of this object and the
   *          DFSphericalother.
   *
   */
  DFSpherical<T> operator/(const DFSpherical<T>& other);
  /** \brief The * operator multiplies the content of this object and the scalar
   *         other and returns the result. 
   *
   *  \return DFSphericalobject which is the product this object and
   *          the the template value T other.
   *
   */
  DFSpherical<T> operator*(const T other);
  /** \brief The * operator for (<typename T>) * DFSpherical. 
   *
   *  \return DFSpherical object which is the product of a the variable
   *          current and the DFSpherical other.
   *
   */
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wnon-template-friend"
  friend DFSpherical<T> operator*(const T current,
      const DFSpherical<T>& other);
  /** \brief Exponential function that returns a DFSpherical where all
   *         elements of DFSpherical other were exponentiated.
   *
   *  \return DFSpherical with exponentiated values of the
   *          DFSpherical other.
   *
   */
  friend DFSpherical<T> exp(const DFSpherical<T>& other);
  /** \brief Natural logarithm that returns a DFSpherical where all
   *         elements of DFSpherical other were logarithmized.
   *
   *  \return DFSpherical with logarithmized values of the
   *          DFSpherical other.
   *
   */
  friend DFSpherical<T> log_natural(const DFSpherical<T>& other);
  #pragma GCC diagnostic pop
  /** \brief std::exception BadSizeException.
   */
  class BadSizeException : public std::exception {
   public:
    /** \brief Overwrite the exception information function what().
     *
     *  \return the text "Size of DFSphericals does not match.".
     */
    virtual const char* what(void) const throw() {
      return "Size of DFSphericals does not match.";
    }
  }
  /** \brief Exception BadSizeException bad_size_error.
   *
   *  The exception is thrown if two DFSphericals do not match in size. 
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
#endif  // SRC_DF_SPHERICAL_HPP_
