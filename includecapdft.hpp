// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_CONVERGENCE_CRITERION_HPP_
#define SRC_CONVERGENCE_CRITERION_HPP_
/** \file src/convergence_criterion.hpp
 *  \brief This file contains the declarations of the ConvergenceCriterion
 *         class.
 */
#include <vector>
#include <string>
#include "data_frame.hpp"  // NOLINT
/** \brief This class defines convergence criteria based on the old and the new
 *         density profiles.
 */
class ConvergenceCriterion {
 public:
  /** \brief Constructor for floating point thresholds
   */
  ConvergenceCriterion(
      const std::vector<DataFrame<1, double>>& old_profile,
      const std::vector<DataFrame<1, double>>& new_profile,
      const double& threshold);
  /** \brief Constructor for integer thresholds
   *
   *  This constructor is used, when counting steps etc.
   */
  ConvergenceCriterion(
      const std::vector<DataFrame<1, double>>& old_profile,
      const std::vector<DataFrame<1, double>>& new_profile,
      const int& threshold_int);
  /** \brief Destructor
   *
   */
  virtual ~ConvergenceCriterion() {}
  /** \brief Check if the criterion is fulfilled (fulfilled = true)
   *
   *  The implementation of the criterion should go here.
   *
   *  \param Uses a pointer to return the convergence progress
   *         (e.g. the largest deviation between the profiles).
   *
   *  \param The threshold determines, when the convergence criterion is reached
   *         given the old and the new density profiles.
   *
   *  \return true if converged, false otherwise
   */
  virtual bool check(double* progress) = 0;
  /** \brief Return name of the criterion
   */
  virtual std::string name() = 0;
  /** \brief Pointer to the old density profile
   */
  const std::vector<DataFrame<1, double>>& old_profile;
  /** \brief Pointer to the proposed density profile
   */
  const std::vector<DataFrame<1, double>>& new_profile;
  /** \brief Threshold value used in the check function (floating point)
   */
  double threshold;
  /** \brief Threshold value used in the check function (integer)
   */
  int threshold_int;
};
#endif  // SRC_CONVERGENCE_CRITERION_HPP_
// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_CONVERGENCE_CRITERION_MAX_DEV_HPP_
#define SRC_CONVERGENCE_CRITERION_MAX_DEV_HPP_
/** \file src/convergence_criterion_max_dev.hpp
 *  \brief This file contains the declarations of the ConvergenceCriterionMaxDev
 *         class.
 */
#include "convergence_criterion.hpp"  // NOLINT
#include <vector>
#include <string>
#include "data_frame.hpp"  // NOLINT
/** \brief This class defines a convergence criterion based on the largest
 *         occuring deviation between old and new density profiles.
 */
class ConvergenceCriterionMaxDev : public ConvergenceCriterion {
 public:
  /** \brief Constructor
   */
  ConvergenceCriterionMaxDev(
      const std::vector<DataFrame<1, double>>& old_profile,
      const std::vector<DataFrame<1, double>>& new_profile,
      const double& target_deviation);
  /** \brief Destructor
   */
  virtual ~ConvergenceCriterionMaxDev();
  /** \brief Check if the criterion is fulfilled (fulfilled = true)
   *
   *  Calculates the largest difference between old and new (proposed)
   *  density profiles and checks if it is smaller than the target_deviation.
   *
   *  \param Uses a pointer to return the current_deviation.
   */
  virtual bool check(double* current_deviation);
  /** \brief Return name of the criterion
   */
  virtual std::string name();
  // The following section is found in the base class
  // /** \brief Pointer to the old density profile
  //  */
  // const std::vector<DataFrame<1, double>>& old_profile;
  // /** \brief Pointer to the proposed density profile
  //  */
  // const std::vector<DataFrame<1, double>>& new_profile;
  // /** \brief Threshold value used in the check function
  //  */
  // double threshold;
  // /** \brief Threshold value used in the check function (integer)
  //  */
  // int threshold_int;
};
#endif  // SRC_CONVERGENCE_CRITERION_MAX_DEV_HPP_
// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_CONVERGENCE_CRITERION_NAN_HPP_
#define SRC_CONVERGENCE_CRITERION_NAN_HPP_
/** \file src/convergence_criterion_nan.hpp
 *  \brief This file contains the declarations of the ConvergenceCriterionNan
 *         class.
 */
#include "convergence_criterion.hpp"  // NOLINT
#include <vector>
#include <string>
#include "data_frame.hpp"  // NOLINT
/** \brief This class defines a convergence criterion based on the number of
 *         nans in the density profiles.
 *
 *  This is usually used to interrupt the iteration process if the iterator
 *  "exploded".
 */
class ConvergenceCriterionNan : public ConvergenceCriterion {
 public:
  /** \brief Constructor
   */
  ConvergenceCriterionNan(
      const std::vector<DataFrame<1, double>>& old_profile,
      const std::vector<DataFrame<1, double>>& new_profile,
      const int& nan_threshold);
  /** \brief Destructor
   */
  virtual ~ConvergenceCriterionNan();
  /** \brief Check if the criterion is fulfilled (fulfilled = true)
   *
   *  Checks for nans.
   *
   *  \param Uses a pointer to return the nan_count.
   */
  virtual bool check(double* current_deviation);
  /** \brief Return name of the criterion
   */
  virtual std::string name();
  // The following section is found in the base class
  // /** \brief Pointer to the old density profile
  //  */
  // const std::vector<DataFrame<1, double>>& old_profile;
  // /** \brief Pointer to the proposed density profile
  //  */
  // const std::vector<DataFrame<1, double>>& new_profile;
  // /** \brief Threshold value used in the check function
  //  */
  // double threshold;
  // /** \brief Threshold value used in the check function (integer)
  //  */
  // int threshold_int;
};
#endif  // SRC_CONVERGENCE_CRITERION_NAN_HPP_
// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_CONVERGENCE_CRITERION_STEPS_HPP_
#define SRC_CONVERGENCE_CRITERION_STEPS_HPP_
/** \file src/convergence_criterion_steps.hpp
 *  \brief This file contains the declarations of the ConvergenceCriterionSteps
 *         class.
 */
#include "convergence_criterion.hpp"  // NOLINT
#include <vector>
#include <string>
#include "data_frame.hpp"  // NOLINT
/** \brief This class defines a convergence criterion based on the number of
 *         iteration steps (i.e. number of times check() was called).
 */
class ConvergenceCriterionSteps : public ConvergenceCriterion {
 public:
  /** \brief Constructor
   */
  ConvergenceCriterionSteps(
      const std::vector<DataFrame<1, double>>& old_profile,
      const std::vector<DataFrame<1, double>>& new_profile,
      const int& steps);
  /** \brief Destructor
   */
  virtual ~ConvergenceCriterionSteps();
  /** \brief Check if the step criterion is fulfilled (fulfilled = true)
   *
   *  Keeps track of the number of steps (number of time check() was called).
   *
   *  \param Uses a pointer to return the number of remaining steps
   */
  virtual bool check(double* current_deviation);
  /** \brief Return name of the criterion
   */
  virtual std::string name();
  // The following section is found in the base class
  // /** \brief Pointer to the old density profile
  //  */
  // const std::vector<DataFrame<1, double>>& old_profile;
  // /** \brief Pointer to the proposed density profile
  //  */
  // const std::vector<DataFrame<1, double>>& new_profile;
  // /** \brief Threshold value used in the check function
  //  */
  // double threshold;
  // /** \brief Threshold value used in the check function (integer)
  //  */
  // int threshold_int;

 private:
  /** \brief Number of steps, i.e. number of times check()
   *
   */
  int step_count;
};
#endif  // SRC_CONVERGENCE_CRITERION_STEPS_HPP_
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
  /** \brief Function that returns the data array pointer.
   *
   */
  T* array();
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
// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_FUNCTIONAL_ES_DELTA_SPHERICAL_HPP_
#define SRC_FUNCTIONAL_ES_DELTA_SPHERICAL_HPP_
/** \file functional_es_delta_spherical.hpp
 *  \brief Header file for the FunctionalESDeltaSpherical class.
 *
 *  The file contains the declarations of the FunctionalESDeltaSpherical class.
 */
// Includes
#include "functional.hpp"  // NOLINT
#include <vector>
#include "data_frame.hpp"  // NOLINT
#include "properties.hpp"  // NOLINT
#include "radial_poisson_solver.hpp"  // NOLINT
/** \brief This class calculates the electrostatic interactions via the delta
 *         functional.
 *  
 *  The theory for this functional was published in
 *  [M. Bültmann and A. Härtel 2022 J. Phys.: Condens. Matter 34 235101].
 */
class FunctionalESDeltaSpherical : public Functional {
 public:
  /** \brief Standard Constructor
   */
  FunctionalESDeltaSpherical();
  /** \brief Manual Constructor
   */
  FunctionalESDeltaSpherical(
      std::vector<DataFrame<1, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties,
      std::vector<size_t> affected_species);
  /** \brief Automated Constructor
   */
  FunctionalESDeltaSpherical(
      std::vector<DataFrame<1, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties);
  /** \brief Destructor
   */
  ~FunctionalESDeltaSpherical();
  /** \brief Calculate the functional derivatives
   *
   *  This function calculates the functional derivatives and updates the
   *  corresponding internal arrays.
   *
   *  It uses the functions:
   *  calc_weighted_densities(),
   */
  virtual void calc_derivative(
      std::vector<DataFrame<1, double>>* functional_derivative);
  /** \brief Calculate bulk derivatives
   */
  virtual void calc_bulk_derivative(std::vector<double>* bulk_derivative);
  /** \brief Calculate the energy value of this functional
   *
   *  Calculate the excess free energy of the delta electrostatic
   *  Functional. For that, one only multiplies the electrostatic potentials of
   *  the weighted densities with the charge densities and integrates the result
   *  over the entire space.
   *
   *  \return Returns the functional energy value
   */
  virtual double calc_energy();

 private:
  /** \brief System length (radius of the sperical geometry)
   */
  double length;
  /** \brief Number of grid points
   */
  size_t grid_count;
  /** \brief Bin width
   */
  double dr;
  /** \brief Bin width in Fourier space
   */
  double dkr;
  /** \brief Number of species
   */
  size_t species_count;
  /** \brief Bjerrum length
   */
  double bjerrum;
  /** \brief System temperature
   */
  double temperature;
  /** \brief Dielectric constant
   */
  double dielectric;
  /** \brief A list with indices of the affected species
   */
  std::vector<size_t> affected_species;
  /** \brief Valency of every species
   *
   *  To be more general, valencies must be given as "double".
   */
  std::vector<double> valencies;
  /** \brief Hard-sphere diameters
   */
  std::vector<double> diameters;
  /** \brief Bulk densities
   */
  std::vector<double> bulk_densities;
  /** \brief Pointer to density profiles
   */
  std::vector<DataFrame<1, double>>* density_profiles_pointer;
  /** \brief Charge density profile for every species
   */
  std::vector<DataFrame<1, double>> charge_density_profiles;
  /** \brief Right hand side of the Poisson equation for each species
   */
  std::vector<DataFrame<1, double>> poisson_rhs;
  /** \brief Electrostatic potentials (numerical solution of Poisson equation)
   *         for each species.
   */
  std::vector<DataFrame<1, double>> potentials;
  /** \brief Fourier transformed delta weight functions
   */
  std::vector<std::vector<DataFrame<1, double>>> weights_delta;
  /** \brief Weighted densities (delta-functions)
   */
  std::vector<std::vector<DataFrame<1, double>>> weighted_densities;
  /** \brief Poisson solver
   *
   *  This object contains the matrix representation of the numerical Poisson
   *  equation.
   */
  RadialPoissonSolver* poisson_solver;
  /** \brief Functional derivatives
   */
  double** functional_derivative;
  /** \brief Extract the system Properties required for this functional
   */
  void extract_system_properties(const Properties& system_properties);
  /** \brief From two of the three electrical properties, the third on can be
   *         calculated.
   */
  void extract_electrical_properties(const Properties& system_properties);
  /** \brief Extract the species specific Properties required for this
   *         functional
   */
  void extract_species_properties(
      const std::vector<Properties>& species_properties);
  /** \brief Allocate memory for all DataFrame objects
   */
  void initialize_all_data_frames();
  /** \brief Initialize the RadialPoissonSolver object
   */
  void initialize_poisson_solver();
  /** \brief Initialize Fourier transformed weight functions
   */
  void initialize_weights();
  /** \brief Calculate the net charge density profile and the rhs of the Poisson
   *         equation.
   */
  void calc_charge_densities();
  /** \brief Calculate the right-hand side of the Poisson equation for every
   *         species.
   */
  void calc_poisson_rhs();
  /** \brief Calculate the net charge density profile and the rhs of the Poisson
   *         equation.
   *
   *  See equation (32a).
   */
  void calc_weighted_densities();
  /** \brief From the charge densities calculate the electrostatic potential
   *
   *  Unlike the mean-field electrostatic functional, the potential of every
   *  species is calculated separately. Moreover, this functional uses weighted
   *  densities instead of the regular charge densities.
   *  The boundary conditions for the solution of the Poisson equations are:
   *  - inner one equals 0 (Neumann) due to the radial symmetry
   *  - outer one equals net charge divided by radial position (Dirichlet)
   *    due to Gauss' theorem.
   *  Note, that there is no external charge at the center.
   *
   */
  void calc_potentials();
  /** \brief Integration over weighted charge densities to obtain net charge
   */
  std::vector<double> calc_charge_weight_dens();

 protected:
};
#endif  // SRC_FUNCTIONAL_ES_DELTA_SPHERICAL_HPP_
// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_FUNCTIONAL_ES_MF_SPHERICAL_HPP_
#define SRC_FUNCTIONAL_ES_MF_SPHERICAL_HPP_
/** \file functional_es_mf_spherical.hpp
 *  \brief Header file for the FunctionalESMFSpherical class.
 *
 *  The file contains the class declarations of the FunctionalESMFSpherical
 *  class.
 *
 */
// Includes
#include "functional.hpp"  // NOLINT
#include <vector>
#include "data_frame.hpp"  // NOLINT
#include "properties.hpp"  // NOLINT
#include "radial_poisson_solver.hpp"  // NOLINT
/** \brief This class calculates the elctrostatic mean field functional.
 *  
 *  This class contains the tools to calculate functional and functional
 *  derivative values for point charges in the mean field approximation.
 *
 */
class FunctionalESMFSpherical : public Functional {
 public:
  /** \brief Standard Constructor
   */
  FunctionalESMFSpherical();
  /** \brief Manual Constructor
   */
  FunctionalESMFSpherical(
      std::vector<DataFrame<1, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties,
      std::vector<size_t> affected_species);
  /** \brief Automated Constructor
   */
  FunctionalESMFSpherical(
      std::vector<DataFrame<1, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties);
  /** \brief Destructor
   */
  ~FunctionalESMFSpherical();
  /** \brief Calculate the functional derivatives
   *
   *  This function calculates the functional derivatives and updates the
   *  corresponding internal arrays.
   *
   *  It uses the functions:
   *  calc_weighted_densities(),
   */
  virtual void calc_derivative(
      std::vector<DataFrame<1, double>>* functional_derivative);
  /** \brief Calculate bulk derivatives
   *
   */
  virtual void calc_bulk_derivative(std::vector<double>* bulk_derivative);
  /** \brief Calculate the energy value of this functional
   *
   *  Calculate the excess free energy of the mean-field electrostatic
   *  Functional. For that, one only multiplies the electrostatic potential with
   *  the total charge density and integrates the result over the entire space.
   *
   *  \return Returns the functional energy value
   */
  virtual double calc_energy();

 private:
  /** \brief System length (radius of the sperical geometry)
   */
  double length;
  /** \brief Number of grid points
   */
  size_t grid_count;
  /** \brief Bin width
   */
  double dr;
  /** \brief Bin width in Fourier space
   */
  double dkr;
  /** \brief Number of species
   */
  size_t species_count;
  /** \brief Bjerrum length
   */
  double bjerrum;
  /** \brief System temperature
   */
  double temperature;
  /** \brief Dielectric constant
   */
  double dielectric;
  /** \brief A list with indices of the affected species
   */
  std::vector<size_t> affected_species;
  /** \brief Valency of every species
   *
   *  To be more general, valencies must be given as "double".
   */
  std::vector<double> valencies;
  /** \brief Bulk densities
   */
  std::vector<double> bulk_densities;
  /** \brief Pointer to density profiles
   */
  std::vector<DataFrame<1, double>>* density_profiles_pointer;
  /** \brief Total charge density profile
   */
  DataFrame<1, double> charge_density_profile;
  /** \brief Right hand side of the Poisson equation
   */
  DataFrame<1, double> poisson_rhs;
  /** \brief Electrostatic potential (numerical solution of Poisson equation)
   */
  DataFrame<1, double> potential;
  /** \brief Poisson solver
   *
   *  This object contains the matrix representation of the numerical Poisson
   *  equation.
   */
  RadialPoissonSolver* poisson_solver;
  /** \brief Functional derivatives */
  double** functional_derivative;
  /** \brief Extract the system Properties required for this functional */
  void extract_system_properties(const Properties& system_properties);
  /** \brief From two of the three electrical properties, the third on can be
   *         calculated.
   */
  void extract_electrical_properties(const Properties& system_properties);
  /** \brief Extract the species specific Properties required for this
   *         functional
   */
  void extract_species_properties(
      const std::vector<Properties>& species_properties);
  /** \brief Allocate memory for all DataFrame objects */
  void initialize_all_data_frames();
  /** \brief Initialize the RadialPoissonSolver object */
  void initialize_poisson_solver();
  /** \brief Calculate the net charge density profile and the rhs of the Poisson
   *         equation.
   */
  void calc_charge_densities();
  /** \brief From the charge densities calculate the electrostatic potential */
  void calc_potential();
  /** \brief Integration over charge densities to obtain net charge */
  double calc_net_charge();

 protected:
};
#endif  // SRC_FUNCTIONAL_ES_MF_SPHERICAL_HPP_
// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2019 Andreas Härtel <http://andreashaertel.anno1982.de/>  // NOLINT
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_FUNCTIONAL_FMT_PLANAR_HPP_
#define SRC_FUNCTIONAL_FMT_PLANAR_HPP_
/** \file functional_fmt_planar.hpp
 *  \brief Header file for the FunctionalFMTPlanar class.
 *
 *  The file contains the class declarations of the FunctionalFMTPlanar
 *  class.
 *
 */
// _____________________________________________________________________________
// Includes
#include <fftw3.h>
#include <vector>
#include "data_frame.hpp"  // NOLINT
#include "functional.hpp"  // NOLINT
#include "properties.hpp"  // NOLINT
// Class forward declarations
// _____________________________________________________________________________
/** \brief FunctionalFMTPlanar calculates the FMT functional in the planar
 *  geometry
 *  
 *  This class contains the tools to calculate functional and functional
 *  derivative values.
 *  It is derived from the abstract Functional base class.
 *
 */
class FunctionalFMTPlanar {
 public:
  /** \brief Standard Constructor
   *
   */
  FunctionalFMTPlanar();
  /** \brief Manual Constructor
   *
   * This constructor chooses the particle species that are supplied from the
   * affected_species vector.
   *
   */
  FunctionalFMTPlanar(
      const std::vector<DataFrame<1, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties,
      const std::vector<size_t>& affected_species);
  /** \brief Automated Constructor
   *
   * This constructor chooses the particle species that are interacting via this
   * functional by checking if the hard sphere diameter exists.
   *
   */
  FunctionalFMTPlanar(
      const std::vector<DataFrame<1, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties);
  /** \brief Destructor
   *
   */
  ~FunctionalFMTPlanar();
  /** \brief Calculate the functional derivatives
   *
   *  This function calculates the functional derivatives and updates the
   *  supplied DataFrame.
   *
   *  It uses the functions:
   *  calc_weighted_densities(),
   *  check_weighted_densities(),
   *  calc_partial_derivatives() (which uses calc_local_partial_derivatives())
   */
  virtual void calc_derivative(
      std::vector<DataFrame<1, double>>* functional_derivative);
  void calc_derivative_warnings(
      std::vector<DataFrame<1, double>>* functional_derivative);
  /** \brief Calculate bulk derivatives
   *
   */
  virtual void calc_bulk_derivative(std::vector<double>* bulk_derivative);
  /** \brief Calculate the (excess free) energy value of this functional
   *
   *  Calculate the energy value of this functional, which corresponds to the
   *  excess free energy contribution of the hard sphere interactions in
   *  equilibrium.
   *
   *  \return Returns the functional energy value
   */
  virtual double calc_energy();

 private:
  /** \brief System length (radius of the sperical geometry)
   *
   */
  double length;
  /** \brief Number of grid points
   *
   */
  size_t grid_count;
  /** \brief Bin sizes in real and Fourier space
   *
   */
  double dz;
  double dkz;
  /** \brief Number of species
   *
   */
  size_t species_count;
  /** \brief Vector that remembers the species, that are affected by this
   *  functional
   */
  std::vector<size_t> affected_species;
  /** \brief Hard sphere diameters
   *
   */
  std::vector<double> diameters;
  /** \brief Bulk densities
   *
   */
  std::vector<double> bulk_densities;
  /** \brief Pointer to density profiles
   *
   */
  const std::vector<DataFrame<1, double>>* density_profiles_pointer;

 protected:
};
#endif  // SRC_FUNCTIONAL_FMT_PLANAR_HPP_
// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_FUNCTIONAL_FMT_SPHERICAL_HPP_
#define SRC_FUNCTIONAL_FMT_SPHERICAL_HPP_
/** \file functional_fmt_spherical.hpp
 *  \brief Header file for the FunctionalFMTSpherical class.
 *
 *  The file contains the class declarations of the FunctionalFMTSpherical
 *  class.
 *
 */
// _____________________________________________________________________________
// Includes
#include <fftw3.h>
#include <vector>
#include "data_frame.hpp"  // NOLINT
#include "functional.hpp"  // NOLINT
#include "properties.hpp"  // NOLINT
// Class forward declarations
// _____________________________________________________________________________
/** \brief FunctionalFMTSpherical calculates the FMT functional in spherical
 *  geometry
 *  
 *  This class contains the tools to calculate functional and functional
 *  derivative values.
 *  It is derived from the abstract Functional base class.
 *
 */
class FunctionalFMTSpherical : public Functional {
 public:
  /** \brief Standard Constructor
   *
   */
  FunctionalFMTSpherical();
  /** \brief Manual Constructor
   *
   * This constructor chooses the particle species that are supplied from the
   * affected_species vector.
   *
   */
  FunctionalFMTSpherical(
      const std::vector<DataFrame<1, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties,
      const std::vector<size_t>& affected_species);
  /** \brief Automated Constructor
   *
   * This constructor chooses the particle species that are interacting via this
   * functional by checking if the hard sphere diameter exists.
   *
   */
  FunctionalFMTSpherical(
      const std::vector<DataFrame<1, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties);
  /** \brief Destructor
   *
   */
  ~FunctionalFMTSpherical();
  /** \brief Calculate the functional derivatives
   *
   *  This function calculates the functional derivatives and updates the
   *  supplied DataFrame.
   *
   *  It uses the functions:
   *  calc_weighted_densities(),
   *  check_weighted_densities(),
   *  calc_partial_derivatives() (which uses calc_local_partial_derivatives())
   */
  virtual void calc_derivative(
      std::vector<DataFrame<1, double>>* functional_derivative);
  void calc_derivative_warnings(
      std::vector<DataFrame<1, double>>* functional_derivative);
  /** \brief Calculate bulk derivatives
   *
   * 
   *
   */
  virtual void calc_bulk_derivative(std::vector<double>* bulk_derivative);
  /** \brief Calculate the (excess free) energy value of this functional
   *
   *  Calculate the energy value of this functional, which corresponds to the
   *  excess free energy contribution of the hard sphere interactions in
   *  equilibrium.
   *
   *  \return Returns the functional energy value
   */
  virtual double calc_energy();

 private:
  /** \brief System length (radius of the sperical geometry)
   *
   */
  double length;
  /** \brief Number of grid points
   *
   */
  size_t grid_count;
  /** \brief Bin sizes in real and Fourier space
   *
   */
  double dr;
  double dkr;
  /** \brief Normalization factors for the Sine/Cosine transforms
   *
   *  Under normal circumstances one could just use dr/2 and dkr/2 as
   *  normalization factors of sine transforms. However, dkr changes for
   *  cosine transforms. The normalizations can be derived from carefully
   *  reading section 2.5.2 ("Real even/odd DFTs").
   *
   *  Normalization for RODFT00: 2 * (grid_count + 1)
   *  Normalization for REDFT00: 2 * ((grid_count+1) - 1) = 2 * grid_count
   *  (Both need to be multiplied by (pi/2).)
   *
   */
  double norm_sin;
  double norm_cos;
  /** \brief Number of species
   *
   */
  size_t species_count;
  /** \brief Vector that remembers the species, that are affected by this
   *  functional
   */
  std::vector<size_t> affected_species;
  /** \brief Hard sphere diameters
   *
   */
  std::vector<double> diameters;
  /** \brief Bulk densities
   *
   */
  std::vector<double> bulk_densities;
  /** \brief Pointer to density profiles
   *
   */
  const std::vector<DataFrame<1, double>>* density_profiles_pointer;
  /** \brief Density profiles times the radial position
   *  
   *  The sine transform of fftw does not require the point r=0, because it is
   *  zero anyways. The cosine transform, however has a non-zero value there.
   *  For this reason the the array is of length grid_count+1 and a sine
   *  transform starts at the second element.
   *
   */
  std::vector<DataFrame<1, double>> density_profiles_times_r;
  /** \brief Fourier transformed density profile
   *
   *  Since it is the sine transformed density_profiles_times_r it also contains
   *  arrays of length grid_count+1.
   *
   */
  std::vector<DataFrame<1, double>> density_profiles_four;
  /** \brief Weighted densities
   *
   * There are three weighted density types: scalar, vectorial, tensorial.
   * There are four arrays containing the four scalar weighted densities.
   * There are two arrays containing the z-components of the two vectorial
   * weighted densities.
   * There are two arrays containing the first (= second) and third element of
   * the tensorial weighted density.
   *
   */
  std::vector<DataFrame<1, double>> scalar_weighted_dens_real;
  std::vector<DataFrame<1, double>> vector_weighted_dens_real;
  std::vector<DataFrame<1, double>> tensor_weighted_dens_real;
  std::vector<DataFrame<1, double>> scalar_weighted_dens_four;
  std::vector<DataFrame<1, double>> vector_weighted_dens_four;
  std::vector<DataFrame<1, double>> tensor_weighted_dens_four;
  /** \brief Weight functions
   *
   *  In the radially symmetric case we only need a few weight functions.
   *  Every vector element contains all weight function of another species.
   *
   */
  std::vector<std::vector<DataFrame<1, double>>> weights_four;
  /** \brief Partial derivatives of the free energy density w.r.t. the
   *  weighted densities
   *
   * There are three partial derivative types: scalar, vectorial, tensorial.
   * There are four arrays containing the four scalar partial derivatives.
   * There are two arrays containing the z-components of the two vectorial
   * partial derivatives.
   * There are two arrays containing the first (= second) and third element of
   * the tensorial partial derivative.
   *
   */
  std::vector<DataFrame<1, double>> scalar_partial_derivative_real;
  std::vector<DataFrame<1, double>> vector_partial_derivative_real;
  std::vector<DataFrame<1, double>> tensor_partial_derivative_real;
  std::vector<DataFrame<1, double>> scalar_partial_derivative_four;
  std::vector<DataFrame<1, double>> vector_partial_derivative_four;
  std::vector<DataFrame<1, double>> tensor_partial_derivative_four;
  /** \brief Terms of the free energy density w.r.t. the weighted densities
   *
   * They are used as dummy for the convolution of the partial derivatives with
   * the corresponding weight functions.
   *
   */
  std::vector<DataFrame<1, double>> scalar_derivative_terms_four;
  std::vector<DataFrame<1, double>> vector_derivative_terms_four;
  std::vector<DataFrame<1, double>> tensor_derivative_terms_four;
  std::vector<DataFrame<1, double>> scalar_derivative_four;
  std::vector<DataFrame<1, double>> vector_derivative_four;
  std::vector<DataFrame<1, double>> tensor_derivative_four;
  /** \brief Flags for the sine and cosine transforms
   *
   *  The first one preserves the input, while the second one might destroy it.
   *
   */
  static const unsigned int flags_keep = FFTW_MEASURE | FFTW_PRESERVE_INPUT;
  static const unsigned int flags_destroy = FFTW_MEASURE;
  /** \brief Calculate the weighted densities
   *
   */
  void calc_weighted_densities();
  /** \brief Check if unphysical values appear in the weighted densities
   *
   */
  void check_weighted_densities();
  /** \brief Calculate the partial derivatives of the free energy densities
   *
   */
  void calc_partial_derivatives();
  /** \brief Calculate the partial derivatives of the free energy densities at
   *  one position
   *
   */
  void calc_local_partial_derivatives(size_t i);
  /** \brief Calculate the weighted partial derivatives of the free energy
   *  densities
   *
   */
  void calc_weighted_partial_derivatives(
      std::vector<DataFrame<1, double>>* functional_derivative);
  /** \brief Calculate the energy density
   *  
   *  \param The index of the position of which the energy density value is
   *  sought.
   *
   *  \return Returns the functional energy density at the specified position.
   *
   */
  double calc_local_energy_density(size_t position);
  /** \brief From the system object extract the system properties
   *
   */
  void extract_system_properties(const Properties& system_properties);
  /** \brief From the system object extract the species properties
   *
   */
  void extract_species_properties(
      const std::vector<Properties>& species_properties);
  /** \brief Initialize all data frame vectors
   *
   */
  void initialize_all_data_frames();
  /** \brief The density profile is changed outside of this functional, thus
   *  we need to update the internal arrays it before calculating something.
   */
  void update_density_times_r();
  /** \brief Calculate the weight functions in Fourier space
   *  For radially symmetric systems only a few weights are actually needed.
   *  The different weighted densities are obtained by introducing the
   *  prefactors in the convolution.
   *
   */
  void calc_weights();
  /** \brief Radial integration
   *
   */
  double radial_integration(double* data, int n, double delta);
};
#endif  // SRC_FUNCTIONAL_FMT_SPHERICAL_HPP_
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
#include "data_frame.hpp"  // NOLINT
// _____________________________________________________________________________
/** \brief Functional is an abstract class
 *
 *  Functional is an abstract class that defines the interfaces for all
 *  explicit functional implementations. 
 *  Accordingly, every excess free energy functional is a child class of 
 *  Functional and must implement the purely virtual functions. 
 *  The template allows to communicate with a functional without knowing about
 *  its particular implementation. 
 */
class Functional {
 public:
  /** Destructor
   *
   */
  virtual ~Functional() {}
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
  virtual void calc_derivative(
      std::vector<DataFrame<1, double>>* functional_derivative) = 0;
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
// SPDX-FileCopyrightText: 2021 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_ITERATOR_HPP_
#define SRC_ITERATOR_HPP_
/** \file src/iterator.hpp
 *  \brief This file contains the declarations of the Iterator class
 *
 *  The class holds the system and all functionals defined on the system. 
 *  It provides methods to iterate the system's density profiles and to perform
 *  calculations in the framework of DFT. 
 *
 */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <vector>
#include <map>
#include "convergence_criterion.hpp"  // NOLINT
#include "data_frame.hpp"  // NOLINT
#include "functional.hpp" // NOLINT
#include "properties.hpp"  // NOLINT
/** \brief Class that holds the system and all functionals and provides DFT 
 *         methods
 *
 *  The class is a container for the physical system which have to be provided 
 *  during declaration. The system particulalry holds all properties and stores
 *  the density profiles. 
 *  The class further provides methods for calculations in the framework of DFT.
 *
 */
class Iterator {
 public:
  /** \brief Constructor
   * 
   *  \param density_profiles The density profiles of different species given as
   *         DataFrame
   *  \param external_potentials The external potentials for each species
   *
   *  Get the density profiles pointer and the external potentials.
   *  The system data is provided to the functionals and is not needed here.
   *
   */
  // TODO(Andreas): add a functional_ideal object
  Iterator(
      std::vector<DataFrame<1, double>>* density_profiles,
      const std::vector<DataFrame<1, double>>& exp_external_potentials,
      const std::vector<Properties>& species_properties);
  /** \brief Destructor
   *
   */
  ~Iterator();
  /** \brief Add an excess functional
   *
   *  This function adds the specified functional to a list (std::vector).
   *  This works for any capdft functional, because they are derived from the
   *  abstract Functional base class.
   *
   *  For example, the function can be used as follows: 
   *  Iterator my_iterator(my_densities, my_external_potentials);
   *  FunctionalFMTSpherical my_functional(...);
   *  my_iterator.add_excess_functional(my_functional);
   *
   */
  void add_excess_functional(Functional* functional);
  /** \brief Remove a functional according to its index. 
   * 
   *  Removes the specified functional from the DFT framework. 
   *  If the index specifies no functional, nothing will happen. 
   * 
   *  \param index Index of the functional that shall be removed.
   *
   *  \return Indicates whether a functional has been removed or not. 
   *
   */
  void remove_excess_functional(size_t index);
  /** \brief Remove all functionals
   *
   */
  void clear_functionals();
  /** \brief Iterates the system densities according to the Picard iteration
   *         method.
   *
   *  \param Mixing factor of the Picard iteration scheme.
   */
  void run_picard(double mixing);
  /** \brief Iterates the system densities according to the Andersen mixing
   *         method.
             https://doi.org/10.1016/j.fluid.2017.03.023
   *
   *  \param Mixing factor of the Picard iteration scheme.
   *  \param memory is the number of density profiles from the past iterations
   *         the scheme will remember.
   */
  void run_anderson(double mixing, size_t memory);
  /** \brief Calculate the excess free energy.
   *
   *  The grand canonical energy of the system is calculated for its current
   *  state without further iteration. 
   *
   *  \return The grand canonical energy in kT per system volume.
   *
   */
  double calculate_excess_free_energy();
  /** \brief Calculate the grand canonical energy of the system. 
   *
   *  The grand canonical energy of the system is calculated for its current
   *  state without further iteration. 
   *
   *  \return The grand canonical energy in kT per system volume.
   *
   */
  double calculate_gc_energy();
  /** \brief Add a convergence criterion
   *
   *  The standard criterion (ConvergenceCriterionMaxDev) is added when a new
   *  Iterator is created.
   *
   */
  template <typename T>
  void add_convergence_criterion(double threshold) {
    ConvergenceCriterion* criterion = 
        new T(*density_profiles, proposed_densities, threshold);
    convergence_criteria.push_back(criterion);
  }
  template <typename U>
  void add_convergence_criterion(int threshold_int) {
    ConvergenceCriterion* criterion =
        new U(*density_profiles, proposed_densities, threshold_int);
    convergence_criteria.push_back(criterion);
  }
  /** \brief Clear convergence criteria
   *
   *  Clears the convergence criteria. This is especially important when
   *  overwriting the standard criteria. If run_*() is executed without a
   *  ConvergenceCriterion, it will run forever.
   *
   */
  void clear_convergence_criteria();

 private:
  /** \brief Density profiles
   *
   */
  std::vector<DataFrame<1, double>>* density_profiles;
  /** \brief Proposed density profiles by Picard iteration
   *
   */
  std::vector<DataFrame<1, double>> proposed_densities;
  /** \brief External potentials
   *
   */
  const std::vector<DataFrame<1, double>>* exp_external_potentials;
  /** \brief Species properties pointer
   *
   *  Species properties are required for the "bulk density" property.
   *
   */
  const std::vector<Properties>* species_properties;
  /** \brief Container holding all functionals that have benn added
   *
   */
  std::vector<Functional*> excess_functionals;
  /** \brief Container holding all convergence criteria to be checked
   *
   */
  std::vector<ConvergenceCriterion*> convergence_criteria;
  /** \brief Functional derivatives
   *
   */
  std::vector<std::vector<DataFrame<1, double>>> functional_derivatives;
  /** \brief Bulk functional derivatives
   *
   */
  std::vector<std::vector<double>> bulk_derivatives;
  /** \brief Check all convergence criteria and return true if one is fulfilled
   *
   */
  bool check_convergence_criteria();
  /** \brief Calculates all possible scalar product combinations of a profile
   *         history
   *
   *  \param profile_history is a vector of density profiles from which all
   *         scalar product combinations are calculated
   *
   *  \return A matrix that contains all scalar products.
   */
  std::vector<std::vector<double>> calc_scalar_products(
      std::vector<std::vector<DataFrame<1, double>>>& profile_history);
  /** \brief Calculates the linear combination of vectors with scalar_products
   *         that minimizes the 2-norm under the condition that the sum of all
   *         coefficients equals 1.
   *
   *  \param scalar_products is a matrix of all scalar products of the vectors
   *         which are used in the linear combination.
   *
   *  \return Linear combination coefficients
   */
  std::vector<double> shortest_linear_combination(
      std::vector<std::vector<double>>& scalar_products);

 protected:
};
/** \brief Function that defines the root search problem that arises from the
 *         Lagrange multiplier formalism in a way that GSL can process it.
 */
int anderson_f(const gsl_vector* x, void* params, gsl_vector* f);
#endif  // SRC_ITERATOR_HPP_
// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2022 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_PROPERTIES_HPP_
#define SRC_PROPERTIES_HPP_
/** \file properties.hpp
 *  \brief Header file for the Properties class.
 *
 *  The file contains the class declarations of the Properties class.
 *
 */
// _____________________________________________________________________________
// Includes
#include <stdexcept>
#include <iostream>
#include <string>
#include <unordered_map>
#include <typeinfo>
// Class forward declarations
// _____________________________________________________________________________
/** \brief Properties class is a container for species or system properties
 *
 */
class Properties {
  /** \brief The System class has full access on the properties. 
   *
   *  In general, the properties must not been changed once they are defined to 
   *  avoid manipulation of calculations during runtime. However, the System 
   *  needs to have access on properties in order to, for instance, update the 
   *  valancies of all species. For this reason, System is declared as a friend.
   *
   */
  // TODO(Andreas): Keine gute idee, wenn System ein template ist.
  //                Warum nicht die update-Funktion public machen?
  // friend class System;
 public:
  /** \brief Constructors
   *
   */
  Properties();
  Properties(const Properties& other);
  /** \brief Destructor
   *
   */
  ~Properties();
  /** \brief Removes all properties
   *
   */
  void clear();
  /** \brief Returns true if property is contained
   *
   */
  bool contains_property(const std::string& property_name) const;
  /** \brief Add a property with an arbitrary data type
   *
   */
  template<typename T>
  void add_property(
      const std::string& property_name, T property_value) {
    TemplateData<T>* new_property = new TemplateData<T>(property_value);
    properties[property_name] = new_property;
  }
  /** \brief Returns a property with an arbitrary data type
   *
   */
  template<typename T>
  bool get_property(const std::string& property_name, T* property_value) const {
    if (contains_property(property_name)) {
      if (typeid(T) != *(properties.at(property_name)->type)) {
        std::cerr << "Properties::get_property(): \"";
        std::cerr << "Error: Requested type and property type not the same.\"";
        std::cerr << std::endl;
        exit(1);
      }
      *property_value = (dynamic_cast<TemplateData<T>*>(
          properties.at(property_name)))->value;
      return true;
    } else {
      throw &missing_property_error;
      return false;
    }
  }
  /** \brief Data class is a universal type class. 
   *
   *  The Data class is used to derive a TemplateData class that stores data of
   *  a certain type (defined by the user via the template declaration). 
   *  For each TemplateData class is also a Data class, different realizations 
   *  (of different data type) of TemplateData can be stored in a container that
   *  holds Data objects. 
   *
   */
  class Data {
   public:
    virtual ~Data() {}
    const std::type_info *type;
  };
  /** \brief TemplateData class is a container to hold an object of arbitrary 
   *         type. 
   *
   *  TemplateData is a container to hold an object of arbitrary type. In order
   *  to allow to store several TemplateData realizations of different type in 
   *  one container (like a vector or map), TemplateData is a child of the Data
   *  class. 
   *
   */
  template <typename T> class TemplateData : public Data {
   public:
    /** \brief Constructor
     *
     */
    explicit TemplateData(T value) {
      type = &(typeid(T));
      this->value = value;
    }
    /** \brief Stores the value of template type
     *
     */
    T value;
  };
  /** \brief Contains all properties.
   *
   */
  std::unordered_map<std::string, Data*> properties;
  /** \brief std::exception MissingPropertyException.
   */
  class MissingPropertyException : public std::exception {
   public:
    /** \brief Overwrite the exception information function what().
     *
     *  \return the text "The adressed property does not exist.".
     */
    virtual const char* what(void) const throw() {
      return "The requested property does not exist.";
    }
  }
  /** \brief Exception MissingPropertyException missing_property_error.
   *
   *  The exception is thrown if a requested property is missing. 
   */
  missing_property_error;

 protected:
  /** \brief Update a property with an arbitrary data type. 
   * 
   *  If the property already exists, its value is updated. 
   *  Otherwise, the property is added.
   * 
   *  \param property_name Name of the property. 
   *  \param property_value Value of the property. 
   *
   *  Example for calling the function: <br>
   *  update_property<double>("my property", 4.5);
   *
   */
  template<typename T>
  void update_property(
      const std::string& property_name, T property_value) {
    auto search = properties.find(property_name);
    if (search != properties.end()) {
      // Property exists: update value
      delete(properties[property_name]);
      TemplateData<T>* new_property = new TemplateData<T>(property_value);
      properties[property_name] = new_property;
    } else {
      // Property does not exist: add value
      this->add_property<T>(property_name, property_value);
    }
  }
};
#endif  // SRC_PROPERTIES_HPP_
// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_RADIAL_POISSON_SOLVER_HPP_
#define SRC_RADIAL_POISSON_SOLVER_HPP_
/** \file radial_poisson_solver.hpp
 *  \brief Header file for the RadialPoissonSolver class.
 *
 *  The file contains the declarations of the RadialPoissonSolver class.
 */
// Includes
#include <cstddef>
// Flag definition for boundary conditions
enum BoundaryFlag {
  NEUMANN_NEUMANN = 0,
  DIRICHLET_DIRICHLET = 1,
  NEUMANN_DIRICHLET = 2,
  DIRICHLET_NEUMANN = 3
};
/** \brief This class contians tools to solve the radial poisson equation
 * 
 *  The numerical Poisson equation can be rewritten into a matrix equation
 *  containing a tridiagonal matrix. These matrices are called "sparse",
 *  because they mostly contain zeros.
 *  The only non-zero elements are: diagonal elements A(i,i), upper diagonal
 *  elements A(i,i+1), lower diagonal elements A(i,i-1).
 *  Hence this matrix can be solved easily and stored in three arrays.
 *
 */
class RadialPoissonSolver {
 public:
  /** \brief Empty Constructor
   *
   */
  RadialPoissonSolver();
  /** \brief Proper Constructor
   * 
   * Passes the size "dim" of the array to be solved, the grid spacing "dr" of
   * the array and the left boundary condition position "inner_distance".
   *
   */
  RadialPoissonSolver(size_t dim, double dr, double inner_distance);
  /** \brief Destructor
   *
   */
  ~RadialPoissonSolver();
  /** \brief Set the radial Laplace matrix with the chosen boundary conditions
   *
   */
  void set_radial_laplacian(BoundaryFlag flag);
  /** \brief Solve the linear equation system: the boundary must be specified
   *  
   *  For example the DIRICHLET_NEUMANN case expects a left boundary value
   *  and a right first derivative value of the sought function.
   *  The DIRICHLET_DIRICHLET case expects a left and right boundary value of
   *  the sought function.
   *
   */
  void solve(double leftBoundaryValue, double rightBoundaryValue, double* rhs,
      double* solution);
  /** \brief Returns the size of the tridiagonal matrix
   *
   */
  size_t size();

 private:
  /** \brief Number of rows/columns of the square tridiagonal matrix
   *
   */
  size_t dim;
  /** Radial bin size (not needed in carthesian case)
   *
   */
  double dr;
  double drdr;
  /** \brief The innermost radial position
   *
   */
  double inner_distance;
  /** \brief Holds the specified boundary conditions
   *
   */
  BoundaryFlag flag;
  /** \brief Boundary condition dependent positions of the first/last bin
   * 
   *  It is either (dr or dr/2) and determined automatically
   *
   */
  double shiftLeft, shiftRight;
  /** \brief Upper off-diagonal elements of the tridiagonal matrix
   *
   */
  double* upper;
  /** \brief Diagonal elements of the tridiagonal matrix
   *
   */
  double* diag;
  /** \brief Lower off-diagonal elements of the tridiagonal matrix
   *
   */
  double* lower;
  /** \brief Set the Laplace matrix with certain boundary conditions
   *
   *  For NN the boundaries should be half a grid point away from the start/end
   *  point (x_{Wall}=x_{-0.5}|x_0|x_1|...|x_N|x_{N+0.5}=x_{Wall}).
   *  For DD the boundaries should be a full grid point away from the start/end
   *  point (x_{Wall}=x_{-1}|x_0|x_1|...|x_N|x_{N+1}=x_Wall).
   *
   */
  void set_radial_laplacian_NN();
  void set_radial_laplacian_DD();
  void set_radial_laplacian_ND();
  void set_radial_laplacian_DN();
  /** \brief Set the Laplace matrix without boundary conditions
   *
   */
  void set_radial_laplacian();
  /** \brief Add the boundary values to the right-hand side
   *
   */
  void set_boundary_values_laplace_radial(double leftBoundaryValue,
      double rightBoundaryValue, double* rhs);
};
#endif  // SRC_RADIAL_POISSON_SOLVER_HPP_
