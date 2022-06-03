// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2021 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_FUNCTIONAL_HPP_
#define SRC_FUNCTIONAL_HPP_
/** \file functional.hpp
 *  \brief Header file for the Functional class.
 *
 *  The file contains the class declarations of the abstract class Functional.
 *
 */
#include <vector>
#include "data_field.hpp"  // NOLINT
// _____________________________________________________________________________
/** \brief Functional is a template class
 *
 *  Functional is an abstract template class that defines the interfaces for all
 *  explicit functional implementations. 
 *  Accordingly, every excess free energy functional is a child class of 
 *  Functional and must implement the purely virtual functions. 
 *  The template allows to communicate with a functional without knowing about
 *  its particular implementation. 
 */
class Functional {
 public:
  virtual ~Functional() = 0;
  /** \brief Calculate the functional derivatives
   *
   *  This function calculates the functional derivatives and stores the result
   *  in the double DataField functional_derivative. 
   *  The DataField must have the correct dimension, otherwise an exception
   *  is thrown. 
   *
   *  \param functional_derivative The double DataField in which the functional
   *         derivative is stored. All contents are overwritten. The dimension 
   *         of the DataField must correspond to the respective functional 
   *         implementation. 
   */
  virtual void calc_derivative(DataField<double>* functional_derivative) = 0;
  /** \brief Calculate bulk derivatives
   *
   *  The function calculates the bulk derivatives and stores results in the 
   *  double vector bulk_derivative. In bulk densities and the functional 
   *  derivatives are described by simple scalar numbers, because they are 
   *  homogeneous in space. 
   *
   *  \param functional_derivative The double vector in which the functional 
   *         bulk derivative is stored. All contents are overwritten. The 
   *         dimension of the vector must correspond to the number of species
   *         set for the functional. 
   */
  virtual void calc_bulk_derivative(std::vector<double>* bulk_derivative) = 0;
  /** \brief Calculate the energy value of this functional
   *
   *  Calculate the energy value of this functional, which approaches the excess
   *  free energy contribution of the hard sphere interactions.
   *
   *  \return Returns the functional energy value. 
   */
  virtual double calc_energy(void) = 0;

 private:
  // This section should stay empty, because this is the template (abstract)
  // functional class.

 protected:
  // This section should stay empty, because this is the template (abstract)
  // functional class.
};
#endif  // SRC_FUNCTIONAL_HPP_
