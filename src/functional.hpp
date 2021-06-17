// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2021 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_FUNCTIONAL_HPP_
#define SRC_FUNCTIONAL_HPP_
/** \file functional.hpp
 *  \brief Header file for the Functional class.
 *
 *  The file contains the class declarations of the Functional 
 *  abstract class.
 *
 */
// _____________________________________________________________________________
/** \brief Functional is a template class
 *
 *  Functional is an abstract template class that dictates the interfaces of all
 *  other functionals. Every excess free energy functional must be a 
 *  child class of Functional and override the purely virtual functions. 
 */
class Functional {
 public:
  /** \brief Constructor
   *
   */
  Functional();
  /** \brief Destructor
   *
   */
  virtual ~Functional();
  /** \brief Calculate the functional derivatives
   *
   *  This function calculates the functional derivatives and updates the
   *  corresponding internal array.
   */
  virtual void calc_derivative() =0;
  /** Calculate bulk derivatives
   *
   */
  virtual void calc_bulk_derivative() =0;
  /** \brief Calculate the energy value of this functional
   *
   *  Calculate the energy value of this functional, which approaches the excess
   *  free energy contribution of the hard sphere interactions.
   *
   *  \return Returns the functional energy value
   */
  virtual double calc_energy() =0;
  /** \brief Calculate the ideal energy (not part of this functional)
   *
   *  Calculates the ideal free energy of all species that interact via this
   *  functional. The ideal free energy ist not really part of this excess
   *  functional.
   *
   *  \return Returns the ideal free energy value
   */
  virtual double calc_ideal_free_energy();  // TODO(Moritz): id gas functional
  /** \brief Calculate the grand potential (not part of this functional)
   *
   *  Calculates the grand potential given a certain chemical potential and
   *  external potential.
   *
   *  \return Returns the grand potential energy
   */
  virtual double calc_grand_potential();  // TODO(Moritz): ideal gas functional

 private:
  // This section should stay empty, because this is the template
  // functional class.

 protected:
  // This section should stay empty, because this is the template
  // functional class.
};
#endif  // SRC_FUNCTIONAL_HPP_
