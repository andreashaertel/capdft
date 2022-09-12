// SPDX-FileCopyrightText: 2021 Philipp Pelagejcev
// <philipp.pelagejcev@physik.uni-freiburg.de>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef SRC_FUNCTIONAL_FMT_CARTESIAN_HPP_
#define SRC_FUNCTIONAL_FMT_CARTESIAN_HPP_

/** \file functional_fmt_cartesian.hpp
 *  \brief Header file of the FunctionalFMTCartesian class.
 *
 *  This header file contains the class declarations of the
 *  FunctionalFMTCartesian class. The FunctionalFMTCartesian class is a
 *  functional that uses fundamental measure theory (FMT) to describe the
 *  interactions of hard spheres (of different sizes) in one dimension. The
 *  functional depends only on the volume fraction. For more information see
 *  White Bear Mark II "Hansen-Goos, H., & Roth, R. (2006) Journal of Physics:
 *  Condensed Matter, 18(37), 8413"
 *
 */

// _____________________________________________________________________________
// Includes
#include <vector>
#include "properties.hpp"
#include "functional.hpp"
#include <fftw3.hpp>
// ____________________________________________________________________________

/** \brief FunctionalFMTCartesian calculates the FMT functional in cartesian
 * geometry
 *
 * This class calculates the FMT functional that describes hard spheres
 * interactions in cartesian geometry. It also contains functions to calculate
 * the functional derivative.
 *
 */

class FunctionalFMTCartesian : public Functional {
 public:
   /** \brief Constructor
    *
    */
  FunctionalFMTCartesian(const Properties& system_properties,
      const std::vector<Properties>& species_properties,
      double** density profile);

  /** \brief Destructor
   *
   */
  ~FunctionalFMTCartesian(void);

  /** \brief Calculates the functional derivative with respect to the density
   *
   * This function calculates the functional derivate with respect to the
   * density of the system. The working steps are split into the following 
   * functions:
   *
   * calc_weighted_densities();
   * calc_partial_derivatives();
   */
  virtual void calc_derivative(void);

  /** \brief Calculate the bulk derivative
   *
   * This function calculates the functional derivative with respect to the
   * density of the system in the bulk limit.
   *
   * (Formulas?)
   *
   */
  virtual void calc_fluid_derivative(void);

  /** \brief This function calculates and returns the energy constribution of this
   * functional to the excess free energy
   *
   * Calculates the energy contribution to the excess free energy and returns
   * the value in dimensions of kT
   *
   */
  virtual double calc_energy(void);

  double **get_derivative_pointer(void);
  double **get_fluid_derivative_pointer(void);

 private:
  /** \brief System length (Distance between electrodes) in nm
   *
   */
  double length;
  size_t grid_count_real;
  size_t grid_count_fourier;
  size_t species_count;
  std::vector<double> diameters;
  std::vector<double> bulk_densities;
  double dz;
  double dkz;
  double** functional_derivative;

  void calc_weighted_densities(void);
  void calc_partial_derivatives(void);
};

#endif  // SRC_FUNCTIONAL_FMT_CARTESIAN_HPP_
