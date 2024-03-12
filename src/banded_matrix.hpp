// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_BANDED_MATRIX_HPP_
#define SRC_BANDED_MATRIX_HPP_
/** \file banded_matrix.hpp
 *  \brief Header file for the BandedMatrix class.
 *
 *  The file contains the class declarations of the BandedMatrix class.
 */
// Includes
#include <cstddef>
#include <vector>
/** \brief This class represents a banded matrix
 *
 *  This class represents a (square) banded matrix with size "size", upper_bw
 *  upper bands and lower_bw lower bands, and of course one diagonal band.
 *
 *  It contains functionality to solve the matrix equations.
 *  For that purpose it uses the GNU Scientific Library.
 *
 *  To use this class for special purposes, like Laplace equations it would be
 *  advisable to derive another class from this class.
 */
class BandedMatrix {
 public:
  /** \brief Empty constructor
   */
  BandedMatrix();
  /** \brief Constructor
   */
  BandedMatrix(size_t row_count, size_t upper_bw, size_t lower_bw);
  /** \brief Destructor
   */
  ~BandedMatrix();
  /** \brief Set all elements to zero
   */
  void zero();
  /** \brief Get the size of the matrix
   */
  size_t size();
  /** \brief Accessor to the band elements
   *  
   *  Get access to the value of the element "band_index" of the "band"th band.
   *  The band_index is equivalent to the row index, which means that there are
   *  forbidden entries, like the first entries in the lower bands and the last
   *  elements in the upper bands.
   */
  double& at(size_t band, size_t index);
  /** \brief Accessor to the matrix element
   *  
   *  Get access to the value of the matrix element A(row,col).
   *  Setting the values in this way does not make sense.
   *  If a non-band element is accessed, an error is thrown.
   */
  double& element(size_t row, size_t col);
  /** \brief Perform the LU decomposition
   *
   *  After specifying the matrix entries this has to be done only once.
   *  The LU decomposition can then be used to solve the matrix equation for
   *  any right-hand side.
   */
  void lu_decomposition();
  /** \brief Solving algorithm
   *
   *  This algorithm uses the LU decomposition for banded matrices from the
   *  GSL library to solve the matrix equation Ax=b, where A is the
   *  banded matrix, x is the vector to be solved for, and b is a known vector.
   */
  void lu_solve(double* rhs, double* solution);

 private: 
  /** \brief Number of rows in the matrix
   */
  size_t row_count;
  /** \brief Number of upper bands (upper band width)
   */
  size_t upper_bw;
  /** \brief Number of lower bands (lower band width)
   */
  size_t lower_bw;
  /** \brief The bands the matrix consists of as a collection of vectors
   */
  std::vector<std::vector<double>> bands;
  /** \brief Upper triagular matrix of the LU decomposition
   */
  std::vector<std::vector<double>> u_matrix;
  /** \brief Lower triagular matrix of the LU decomposition
   */
  std::vector<std::vector<double>> l_matrix;
  /** \brief Keeps track of the pivoting in the lu_decomposition()
   */
  std::vector<size_t> pivot;
  /** \brief Check if the element to be accessed is a member of one of the bands
   */
  bool is_out_of_bounds(size_t band, size_t index);
  /** \brief Free the memory and recreate the banded matrix
   */
  void recreate_banded_matrix();
  /** \brief Free the memory and recreate all containers needed for the LU
   *         decomposition
   */
  void recreate_lu_containers();
};
#endif  // SRC_BANDED_MATRIX_HPP_
