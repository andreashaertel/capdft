// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_SPARSE_MATRIX_HPP_
#define SRC_SPARSE_MATRIX_HPP_
/** \file sparse_matrix.hpp
 *  \brief Header file for the SparseMatrix class.
 *
 *  The file contains the class declarations of the SparseMatrix class.
 */
// Includes
#include <list>
#include <vector>
/** \brief This class represents a general sparse matrix
 *
 *  This class represents a general sparse matrix.
 *  The matrix is represented by 2 vectors.
 *  row_index keeps track of the row position of the non-zero
 *  values, while values contains the non-zero values of the matrix.
 *  Both vectors have col_count many elements (1 for every column).
 *  The storage format is based on the widely known 
 *  "compressed column storage format".
 *
 *  This class contains functionality to solve the matrix equations iteratively.
 *
 *  To use this class for special purposes, like Laplace equations it would be
 *  advisable to derive another class from this one.
 */
class SparseMatrix {
 public:
  /** \brief Constructor
   */
  SparseMatrix(size_t row_count, size_t col_count, size_t row_reserved_space);
  /** \brief Constructor
   */
  SparseMatrix(size_t row_count, size_t col_count);
  /** \brief Constructor
   */
  SparseMatrix(size_t row_col_count);
  /** \brief Empty constructor
   */
  SparseMatrix();
  /** \brief Destructor
   */
  ~SparseMatrix();
  /** \brief Set all elements to zero
   */
  void zero();
  /** \brief Get the value of a matrix element
   *
   *  Also returns a zero if an element is accessed, that is not explicitely
   *  contained in the values vector.
   */
  double get(size_t row, size_t col);
  /** \brief Set the value of a matrix element
   *  
   *  Also allocates memory if a zero element is changed to a non-zero element.
   */
  void set(size_t row, size_t col, double value);
  /** \brief Multiply this matrix to a vector (A*x)
   */
  std::vector<double> matrix_times_vector(std::vector<double>& x);
  /** \brief Multiply a transposed vector by this matrix (x*A)
   */
  std::vector<double> vector_times_matrix(std::vector<double>& x);
  /** \brief Scalar product for completeness; matrix is not involved
   */
  double vector_times_vector(std::vector<double>& x, std::vector<double>& y);
  /** \brief Find a solution to the matrix equation A*x=b
   *
   *  The solution algorithm is called GMRES.
   *  \return True if the algorithm converged within the specified epsilon.
   */
  bool solve(
      std::vector<double>& rhs, std::vector<double>& solution,
      size_t max_iterations, double epsilon);

 private: 
  /** \brief Number of rows in the matrix
   */
  size_t row_count;
  /** \brief Number of columns in the matrix
   */
  size_t col_count;
  /** \brief The number of non-zero elements in the matrix
   */
  size_t size;
  /** \brief The values of the sparse matrix
   */
  std::vector<std::vector<double>> values;
  /** \brief The row indices of the values of the sparse matrix
   */
  std::vector<std::vector<size_t>> row_index;
  /** \brief Erase an non-zero element
   */
  void erase(size_t position, size_t col);
  /** \brief Insert a non-zero element
   */
  void insert(size_t position, size_t row, size_t col, double value);
};
#endif  // SRC_SPARSE_MATRIX_HPP_
