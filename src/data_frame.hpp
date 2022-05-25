// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2022 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_DATA_FRAME_HPP_
#define SRC_DATA_FRAME_HPP_
/** \file data_frame.hpp
 *  \brief This file contains the declarations of the virtual DataFrame class
 *
 */
#include "properties.hpp"
/** \brief Container class for general data (e.g. density profiles or functional
 *         derivatives)
 *
 *  The DataFrame is the general data object of the capDFT project which can 
 *  hold data profiles as density profiles, functional derivatives, or external
 *  potential fields. The DataFrame provides functionality required by classes 
 *  like System or Dft to perform calculations on data (without any knowledge
 *  on the specific implementation). 
 * 
 *  Specific implementations of a DataFrame have to be derived from DataFrame 
 *  and must agree with the respective functional implementations of this 
 *  specific geometry. A name convention for derived realizations is: <br>
 *  Let the specific geometry be called GeometryA. Then the respective 
 *  DataFrame child should be called DataFrameGeometryA and be stored in the 
 *  files data_frame_geometry_a.*. The respective functional implementations 
 *  should follow a corresponding name convention, functional_KIND_geometry_a.*.
 *
 *  A DataFrame must implement the following operators: 
 *   -# DataFrame = DataFrame   (virtual)
 *   -# DataFrame += DataFrame  (virtual)
 *   -# DataFrame -= DataFrame  (virtual)
 *   -# DataFrame *= DataFrame  (virtual)
 *   -# DataFrame /= DataFrame  (virtual)
 *   -# DataFrame *= (double)   (virtual)
 *   -# DataFrame + DataFrame   (virtual)
 *   -# DataFrame - DataFrame   (virtual)
 *   -# DataFrame * DataFrame   (virtual)
 *   -# DataFrame / DataFrame   (virtual)
 *   -# DataFrame * (double)    (virtual)
 *   -# (double) * DataFrame
 *   -# exp (DataFrame)         (virtual)
 *
 */
class DataFrame {
 public:
  /** \brief Virtual Constructors
   *
   *  A DataFrame must have two Constructors: 
   *  One to construct a new implementation, and a second one to create a copy 
   *  of an existing one. To create a new implementation, parameters may be 
   *  required: They can be provided via a Properties object, for instance the 
   *  system properties of the System class. 
   *  The Properties must must be offered in the Constructor.
   *
   *  The copy Constructor, however, just creates a copy of an existing 
   *  DataFrame realization and, thus, does not require the Template type 
   *  declaration. 
   *
   */
  DataFrame(const Properties properties);
  DataFrame(const DataFrame& other);
  /** \brief Virtual destructor
   *
   */
  virtual ~DataFrame() = 0;
  /** \brief Test for DataFrame other having the same size as this. 
   *
   *  \param other The DataFrame thats size is compared to this DataFrame. 
   *
   *  \return True, if the other DataFrame has the same size as this one, 
   *          False otherwise. 
   *
   */
  virtual bool same_size(const DataFrame& other) = 0;

  /** \brief The = operator copies the content of other into this. */
  virtual DataFrame& operator=(const DataFrame& other) = 0;
  /** \brief The += operator adds the content of other to this and then returns
   *         this. 
   */
  virtual DataFrame& operator+=(const DataFrame& other) = 0;
  /** \brief The -= operator subtracts the content of other from this and then
   *         returns this. 
   */
  virtual DataFrame& operator-=(const DataFrame& other) = 0;
  /** \brief The *= operator multiplies the content of other to this and then
   *         returns this. 
   */
  virtual DataFrame& operator*=(const DataFrame& other) = 0;
  /** \brief The /= operator divides this by the content of other and then 
   *         returns this. 
   */
  virtual DataFrame& operator/=(const DataFrame& other) = 0;
  /** \brief The *= operator for scalar multiplication multiplies all entries of
   *         this with the scalar value other and then returns this. 
   */
  virtual DataFrame& operator*=(const double other) = 0;
  /** \brief The + operator adds the content of this and other and returns the 
   *         result. 
   *
   *  Use the following code: 
   *
   *      DataFrame& operator+(const DataFrame& other) {
   *        // Check for correct sizes
   *        if (! this->same_size(other))
   *        throw &bad_size_error;
   *        // Use existing copy Constructor and += operators.
   *        DataFrame result(*this);
   *        result += other;
   *        return result;
   *      };
   *
   */
  virtual DataFrame& operator+(const DataFrame& other) = 0;
  /** \brief The - operator subtracts the content of this and other and returns
   *         the result. 
   *
   *  Use the following code: 
   *
   *      DataFrame& DataFrame::operator-(const DataFrame& other) {
   *        // Check for correct sizes
   *        if (! this->same_size(other))
   *        throw &bad_size_error;
   *        // Use existing copy Constructor and -= operators.
   *        DataFrame result(*this);
   *        result -= other;
   *        return result;
   *      }
   *
   */
  virtual DataFrame& operator-(const DataFrame& other) = 0;
  /** \brief The * operator multiplies the content of this and other and returns
   *         the result. 
   *
   *  Use the following code: 
   *
   *      DataFrame& DataFrame::operator*(const DataFrame& other) {
   *        // Check for correct sizes
   *        if (! this->same_size(other))
   *        throw &bad_size_error;
   *        // Use existing copy Constructor and *= operators.
   *        DataFrame result(*this);
   *        result *= other;
   *        return result;
   *      }
   *
   */
  virtual DataFrame& operator*(const DataFrame& other) = 0;
  /** \brief The / operator divides the content of this and other and returns 
   *         the result. 
   *
   *  Use the following code: 
   *
   *      DataFrame& DataFrame::operator/(const DataFrame& other) {
   *      // Check for correct sizes
   *      if (! this->same_size(other))
   *      throw &bad_size_error;
   *      // Use existing copy Constructor and /= operators.
   *      DataFrame result(*this);
   *      result /= other;
   *      return result;
   *    }
   *
   */
  virtual DataFrame& operator/(const DataFrame& other) = 0;
  /** \brief The * operator multiplies the content of this and the scalar other
   *         and returns the result. 
   *
   *  Use the following code: 
   *
   *      DataFrame& DataFrame::operator*(const double other) {
   *        // Use existing copy Constructor and *= for double operators.
   *        DataFrame result(*this);
   *        result *= other;
   *        return result;
   *      }
   *
   */
  virtual DataFrame& operator*(const double other) = 0;
  /** \brief The * operator for (double) * DataFrame. 
   *
   *  Use the following code: 
   *
   *      friend DataFrame& operator*(const double current, const DataFrame& 
   *        other) {
   *        // Use existing copy Constructor and *= for double operators.
   *        DataFrame result(other);
   *        result *= current;
   *        return result;
   *      }
   *
   */
  friend DataFrame& operator*(const double current, const DataFrame& other);
  /** \brief Exponential function applied to all elements of DataFrame. 
   *
   */
  virtual DataFrame& exp() = 0;
  /** \brief Natural logarithm function applied to all elements of DataFrame. 
   *
   */
  virtual DataFrame& log_natural() = 0;
  /** \brief std::exception BadSizeException.
   */
  class BadSizeException : public std::exception {
   public:
    /** \brief Overwrite the exception information function what().
     *
     *  \return the text "Size of DataFrame's does not match.".
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

 private: 

};
#endif  // SRC_DATA_FRAME_HPP_
