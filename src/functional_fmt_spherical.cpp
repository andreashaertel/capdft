// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "functional_fmt_spherical.hpp"
#include <fftw3.h>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
// _____________________________________________________________________________
FunctionalFMTSpherical::FunctionalFMTSpherical() {
  //
}
// _____________________________________________________________________________
FunctionalFMTSpherical::FunctionalFMTSpherical(
    System* system, const std::vector<size_t>& affected_species)
  : affected_species(affected_species) {
  // Clear all std::vectors
  diameters.clear();
  bulk_densities.clear();
  // Get system properties
  extract_system_properties(system);
  // Get species properties; excludes all species without diameter property
  extract_species_properties(system);
  // Initialize the density profile and density profile times radial position
  density_profile_pointer = system->get_density_profile_pointer();
  density_profile_times_r = new DataField<double>(species_count, grid_count+1);
  density_profile_four = new DataField<double>(species_count, grid_count+1);
  update_density_times_r();
  // Initialize weighted densities
  scalar_weighted_dens_real = new DataField<double>(4, grid_count);
  vector_weighted_dens_real = new DataField<double>(2, grid_count);
  tensor_weighted_dens_real = new DataField<double>(2, grid_count);
  scalar_weighted_dens_four = new DataField<double>(4, grid_count+1);
  vector_weighted_dens_four = new DataField<double>(4, grid_count+1);
  tensor_weighted_dens_four = new DataField<double>(3, grid_count+1);
  // Initialize partial derivatives
  scalar_partial_derivative_real = new DataField<double>(4, grid_count);
  vector_partial_derivative_real = new DataField<double>(2, grid_count);
  tensor_partial_derivative_real = new DataField<double>(2, grid_count);
  scalar_partial_derivative_four = new DataField<double>(4, grid_count+1);
  vector_partial_derivative_four = new DataField<double>(4, grid_count+1);
  tensor_partial_derivative_four = new DataField<double>(4, grid_count+1);
  // Initialize derivative terms
  scalar_derivative_four = new DataField<double>(4, grid_count);
  vector_derivative_four = new DataField<double>(2, grid_count);
  tensor_derivative_four = new DataField<double>(2, grid_count);
  scalar_derivative_terms_four = new DataField<double>(4, grid_count+1);
  vector_derivative_terms_four = new DataField<double>(4, grid_count+1);
  tensor_derivative_terms_four = new DataField<double>(4, grid_count+1);
  // Initialize weights
  for (size_t i = 0; i != species_count; ++i) {
    weights_four.push_back(DataField<double>(3, grid_count+1));
  }
  initialize_weights();
}
// _____________________________________________________________________________
FunctionalFMTSpherical::FunctionalFMTSpherical(System* system)
  : FunctionalFMTSpherical(system, std::vector<size_t>(0)) {
}
// _____________________________________________________________________________
FunctionalFMTSpherical::~FunctionalFMTSpherical() {
  // Free memory of density profiles
  density_profile_times_r->~DataField();
  density_profile_four->~DataField();
  // Free memory of weighted densities
  scalar_weighted_dens_real->~DataField();
  vector_weighted_dens_real->~DataField();
  tensor_weighted_dens_real->~DataField();
  scalar_weighted_dens_four->~DataField();
  vector_weighted_dens_four->~DataField();
  tensor_weighted_dens_four->~DataField();
  // Free memory of partial derivatives
  scalar_partial_derivative_real->~DataField();
  vector_partial_derivative_real->~DataField();
  tensor_partial_derivative_real->~DataField();
  scalar_partial_derivative_four->~DataField();
  vector_partial_derivative_four->~DataField();
  tensor_partial_derivative_four->~DataField();
  // Free memory of derivative dummy variablea
  scalar_derivative_four->~DataField();
  vector_derivative_four->~DataField();
  tensor_derivative_four->~DataField();
  scalar_derivative_terms_four->~DataField();
  vector_derivative_terms_four->~DataField();
  tensor_derivative_terms_four->~DataField();
}
// _____________________________________________________________________________
void FunctionalFMTSpherical::extract_system_properties(System* sys) {
  Properties properties(sys->get_system_properties());
  properties.get_property("length", &length);
  properties.get_property("grid count", &grid_count);
  // Calculate bin sizes
  dr = length / static_cast<double>(grid_count);
  dkr = (2. * M_PI) / (2. * static_cast<double>(grid_count+1) * dr);
  // Calculate normalization factors for sine and cosine transforms
  norm_sin =
      sqrt(M_PI / (2. * static_cast<double>(2 * (grid_count + 1))));
  norm_cos =
      sqrt(M_PI / (2. * static_cast<double>(2 * grid_count)));
}
// _____________________________________________________________________________
void FunctionalFMTSpherical::extract_species_properties(System* sys) {
  // Sort the affected species numbers
  std::sort(affected_species.begin(), affected_species.end());
  // Remove duplicates
  affected_species.erase(
      unique(affected_species.begin(), affected_species.end()),
      affected_species.end());
  // If no affected species were specified, find them automatically
  std::vector<Properties> spec_prop = sys->get_species_properties();
  double diameter{0.};
  double bulk_density{0.};
  if (affected_species.empty()) {
    for (auto it = spec_prop.begin(); it != spec_prop.end(); ++it) {
      if (it->get_property("diameter", &diameter)) {  // only species with diam.
        if (!it->get_property("bulk density", &bulk_density)) {
          std::cerr << "FunctionalFMTSpherical::extract_species_properties(): ";
          std::cerr << "\"Error: A species with a diameter but no density ";
          std::cerr << "was detected.\"" << std::endl;
          exit(1);
        }
        affected_species.push_back(it - spec_prop.begin());
      }
    }
  }
  // Extract properties
  for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
    if (!spec_prop.at(*it).get_property("diameter", &diameter) ||
        !spec_prop.at(*it).get_property("bulk density", &bulk_density)) {
      std::cerr << "FunctionalFMTSpherical::extract_species_properties(): ";
      std::cerr << "\"Error: One species is missing a required parameter.";
      std::cerr << std::endl;
      exit(1);
    }
    diameters.push_back(diameter);
    bulk_densities.push_back(bulk_density);
  }
  // Count species that interact via the hard sphere potential
  species_count = affected_species.size();
}
// _____________________________________________________________________________
void FunctionalFMTSpherical::update_density_times_r() {
  double r{0.};
  size_t index{0};
  for (size_t i = 0; i < species_count; ++i) {
    index = affected_species.at(i);
    // The following construction is due to the fact, that r*rho(r) contains
    // r = 0 and rho(r) does not.
    for (size_t j = 0; j < grid_count + 1; ++j) {
      r = dr * static_cast<double>(j);
      if (j == 0) {
        density_profile_times_r->at(0, j) = 0.;
      } else {
        density_profile_times_r->at(i, j) =
            r * density_profile_pointer->element(index, j-1);
      }
    }
  }
}
// _____________________________________________________________________________
void FunctionalFMTSpherical::initialize_weights() {
  // Declare auxiliary variables. R means radius of the hard sphere.
  // kr is the radial position in Fourier space. Hence, RR is the square of R
  // and so on.
  double kr{0.}, krkr{0.}, krkrkr{0.};
  double R{0.}, RR{0.}, RRR{0.};
  double Rkr{0.};
  for (size_t i = 0; i != species_count; ++i) {
    R = diameters.at(i) / 2.;
    RR = R * R;
    RRR = RR * R;
    for (size_t j = 0; j != grid_count + 1; ++j) {
      kr = dkr * static_cast<double>(j);
      Rkr = R * kr;
      krkr = kr * kr;
      krkrkr = krkr * kr;
      // Calculate limits kr --> 0
      if (j == 0) {
        // scalar weight 3 (theta)
        weights_four.at(i).at(0, j) = 4. * M_PI * RRR / 3.;
        // scalar weight 2 (delta)
        weights_four.at(i).at(1, j) = 4. * M_PI * RR;
        // tensorial
        weights_four.at(i).at(2, j) =
            weights_four.at(i).at(1, j) - 3. * weights_four.at(i).at(0, j) / R;
      } else {
        // scalar weight 3 (theta)
        weights_four.at(i).at(0, j) = 4. * M_PI * (sin(Rkr) - Rkr * cos(Rkr)) /
            krkrkr;
        // scalar weight 2 (delta)
        weights_four.at(i).at(1, j) = 4. * M_PI * R * sin(Rkr) / kr;
        // tensorial
        weights_four.at(i).at(2, j) =
            weights_four.at(i).at(1, j) - 3. * weights_four.at(i).at(0, j) / R;
      }
    }
  }
}
// _____________________________________________________________________________
void FunctionalFMTSpherical::calc_derivative(
    DataField<double>* functional_derivative) {
  // Check if the given DataField has correct array length
  if (functional_derivative->get_array_size() != grid_count) {
    std::cerr << "FunctionalFMTSpherical::calc_derivative(): \"";
    std::cerr << "Error: Supplied DataField has incorrect array size.\"";
    std::cerr << std::endl;
    exit(1);
  }
  // Calculate the weighted densities
  calc_weighted_densities();
  // Check whether there are unphysical values in the  weighted densities
  check_weighted_densities();
  // From the weighted densities calculate the partial derivatives of the
  // excess free energy
  calc_partial_derivatives();
  // Calculate the weighted partial derivatives
  calc_weighted_partial_derivatives(functional_derivative);
}
// _____________________________________________________________________________
void FunctionalFMTSpherical::calc_derivative_no_warnings(
    DataField<double>* functional_derivative) {
  // Check if the given DataField has correct array length
  if (functional_derivative->get_array_size() != grid_count) {
    std::cerr << "FunctionalFMTSpherical::calc_derivative(): \"";
    std::cerr << "Error: Supplied DataField has incorrect array size.\"";
    std::cerr << std::endl;
    exit(1);
  }
  // Calculate the weighted densities
  calc_weighted_densities();
  // From the weighted densities calculate the partial derivatives of the
  // excess free energy
  calc_partial_derivatives();
  // Calculate the weighted partial derivatives
  calc_weighted_partial_derivatives(functional_derivative);
}
// _____________________________________________________________________________
void FunctionalFMTSpherical::calc_bulk_derivative(
    std::vector<double>* bulk_derivative) {
  // Auxiliary variables
  double diameter{0.}, diameter2{0.}, diameter3{0.};
  // Set all weighted densities to 0
  scalar_weighted_dens_real->zeros();
  // Calculate the fluid weighted densities and save them to the first position
  for (size_t i = 0; i != species_count; ++i) {
    // Define auxiliary values
    diameter = diameters.at(i);
    diameter2 = diameter * diameter;
    diameter3 = diameter2 * diameter;
    // All non-scalar weighted densities are 0 the others are obviously constant
    scalar_weighted_dens_real->at(0, 0) +=
        bulk_densities.at(i) * diameter3 * M_PI / 6.0;
    scalar_weighted_dens_real->at(1, 0) +=
        bulk_densities.at(i) * diameter2 * M_PI;
    scalar_weighted_dens_real->at(2, 0) +=
        bulk_densities.at(i) * diameter / 2.0;
    scalar_weighted_dens_real->at(3, 0) += bulk_densities.at(i);
  }
  // Evaluate the partial derivatives at the first position
  calc_local_partial_derivatives(0);
  std::fill(bulk_derivative->begin(), bulk_derivative->end(), 0.);
  size_t spec_i{0};
  for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
    spec_i = it - affected_species.begin();
    // Define auxiliary values
    diameter = diameters.at(spec_i);
    diameter2 = diameter * diameter;
    diameter3 = diameter2 * diameter;
    // Calculate the derivatives of the excess free energy functional w.r.t. the
    // weighted densities.
    bulk_derivative->at(*it) = scalar_partial_derivative_real->at(0, 0) *
        diameter3 * M_PI / 6.0;
    bulk_derivative->at(*it) += scalar_partial_derivative_real->at(1, 0) *
        diameter2 * M_PI;
    bulk_derivative->at(*it) += scalar_partial_derivative_real->at(2, 0) *
        diameter / 2.0;
    bulk_derivative->at(*it) += scalar_partial_derivative_real->at(3, 0);
  }
}
// _____________________________________________________________________________
double FunctionalFMTSpherical::calc_energy() {
  double integral{0.};
  DataField<double> free_energy_density(1, grid_count);
  // Calculate the weighted densities for the current density profile
  calc_weighted_densities();
  check_weighted_densities();
  // Store free energy density
  for (size_t i = 0; i < grid_count; ++i) {
    free_energy_density.at(0, i) = calc_local_energy_density(i);
  }
  integral = radial_integration(free_energy_density.array(0), grid_count, dr);
  return integral;
}
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
void FunctionalFMTSpherical::calc_weighted_densities() {
  // TODO(Moritz): speed enhancement by not using "at" and by saving transforms
  // TODO(Moritz): speed enhancement by avoiding some transforms
  double r{0.}, rr{0.}, rrr{0.};
  double kr{0.}, krkr{0.};
  std::vector<fftw_plan> forward_plans;
  std::vector<fftw_plan> backward_plans;
  // Specify the plans
  for (size_t i = 0; i != species_count; ++i) {
    forward_plans.push_back(
        fftw_plan_r2r_1d(grid_count,
            density_profile_times_r->array(i) + 1,
            density_profile_four->array(i) + 1,
            FFTW_RODFT00, flags_keep));
  }
  for (size_t i = 0; i != 4; ++i) {  // TODO(Moritz): make number static
    backward_plans.push_back(
        fftw_plan_r2r_1d(grid_count,
            scalar_weighted_dens_four->array(i) + 1,
            scalar_weighted_dens_four->array(i) + 1,
            FFTW_RODFT00, flags_destroy));
  }
  backward_plans.push_back(
      fftw_plan_r2r_1d(grid_count,
          vector_weighted_dens_four->array(0) + 1,
          vector_weighted_dens_four->array(0) + 1,
          FFTW_RODFT00, flags_destroy));
  backward_plans.push_back(
      fftw_plan_r2r_1d(grid_count + 1,
          vector_weighted_dens_four->array(1),
          vector_weighted_dens_four->array(1),
          FFTW_REDFT00, flags_destroy));
  backward_plans.push_back(
      fftw_plan_r2r_1d(grid_count,
          vector_weighted_dens_four->array(2) + 1,
          vector_weighted_dens_four->array(2) + 1,
          FFTW_RODFT00, flags_destroy));
  backward_plans.push_back(
      fftw_plan_r2r_1d(grid_count + 1,
          vector_weighted_dens_four->array(3),
          vector_weighted_dens_four->array(3),
          FFTW_REDFT00, flags_destroy));
  backward_plans.push_back(
      fftw_plan_r2r_1d(grid_count,
          tensor_weighted_dens_four->array(0) + 1,
          tensor_weighted_dens_four->array(0) + 1,
          FFTW_RODFT00, flags_destroy));
  backward_plans.push_back(
      fftw_plan_r2r_1d(grid_count + 1,
          tensor_weighted_dens_four->array(1),
          tensor_weighted_dens_four->array(1),
          FFTW_REDFT00, flags_destroy));
  backward_plans.push_back(
      fftw_plan_r2r_1d(grid_count,
          tensor_weighted_dens_four->array(2) + 1,
          tensor_weighted_dens_four->array(2) + 1,
          FFTW_RODFT00, flags_destroy));
  // Creating plans clears the input array hence we need to update it
  update_density_times_r();
  // Sine transform of r*rho(r)
  for (auto it = forward_plans.begin(); it != forward_plans.end(); ++it) {
    fftw_execute(*it);
  }
  // Make sure that the point r=0 is zero
  for (size_t i = 0; i != species_count; ++i) {
    density_profile_four->at(i, 0) = 0.;
  }
  // Normalize
  for (size_t i = 1; i != grid_count+1; ++i) {
    kr = dkr * static_cast<double>(i);
    for (size_t j = 0; j != species_count; ++j) {
      density_profile_four->at(j, i) *= 4. * M_PI / kr;
      density_profile_four->at(j, i) *= norm_sin;
    }
  }
  // Convolution of density profiles with weights
  scalar_weighted_dens_four->zeros();
  vector_weighted_dens_four->zeros();
  tensor_weighted_dens_four->zeros();
  for (size_t i = 1; i != grid_count+1; ++i) {
    kr = dkr * static_cast<double>(i);
    krkr = kr * kr;
    for (size_t j = 0; j != species_count; ++j) {
      // Scalar weighted density n_3
      scalar_weighted_dens_four->at(0, i) += density_profile_four->at(j, i) *
          kr * weights_four.at(j).at(0, i) * norm_sin;
      // Scalar weighted density n_2
      scalar_weighted_dens_four->at(1, i) += density_profile_four->at(j, i) *
          kr * weights_four.at(j).at(1, i) * norm_sin;
      // Scalar weighted density n_1
      scalar_weighted_dens_four->at(2, i) += density_profile_four->at(j, i) *
          kr * weights_four.at(j).at(1, i) / (2. * M_PI * diameters.at(j)) *
          norm_sin;
      // Scalar weighted density n_0
      scalar_weighted_dens_four->at(3, i) += density_profile_four->at(j, i) *
          kr * weights_four.at(j).at(1, i) / (M_PI * pow(diameters.at(j), 2)) *
          norm_sin;
      // Vectorial weighted density n_2 z-component in rotated system
      // The radial FT of the weight w2 transforms it into the scalar FT w3
      vector_weighted_dens_four->at(0, i) += density_profile_four->at(j, i) *
          kr * weights_four.at(j).at(0, i) * norm_sin;  // sine transform
      vector_weighted_dens_four->at(1, i) += density_profile_four->at(j, i) *
          krkr * weights_four.at(j).at(0, i) * norm_cos;  // cosine transform
      // Vectorial weighted density n_1 z-component in rotated system
      vector_weighted_dens_four->at(2, i) += density_profile_four->at(j, i) *
          kr * weights_four.at(j).at(0, i) * norm_sin /
          (2. * M_PI * diameters.at(j));  // sine transform
      vector_weighted_dens_four->at(3, i) += density_profile_four->at(j, i) *
          krkr * weights_four.at(j).at(0, i) * norm_cos /
          (2. * M_PI * diameters.at(j));  // cosine transform
      // Tensorial weighted density n_m2: the diagonal elements only contain
      // these three terms; all off-diagonal elements are zero
      tensor_weighted_dens_four->at(0, i) += density_profile_four->at(j, i) *
          weights_four.at(j).at(2, i) * norm_sin / kr;  // sin/kr
      tensor_weighted_dens_four->at(1, i) += density_profile_four->at(j, i) *
          weights_four.at(j).at(2, i) * norm_cos;  // cos
      tensor_weighted_dens_four->at(2, i) += density_profile_four->at(j, i) *
          weights_four.at(j).at(2, i) * kr * norm_sin;  // sin*kr
    }
  }
  // Sine and cosine transforms
  for (auto it = backward_plans.begin(); it != backward_plans.end(); ++it) {
    fftw_execute(*it);
  }
  // Normalize
  *(scalar_weighted_dens_four) *= 4. * M_PI;
  *(vector_weighted_dens_four) *= 4. * M_PI;
  *(tensor_weighted_dens_four) *= 4. * M_PI;
  *(scalar_weighted_dens_four) *= 1. / pow(2. * M_PI, 3);
  *(vector_weighted_dens_four) *= 1. / pow(2. * M_PI, 3);
  *(tensor_weighted_dens_four) *= 1. / pow(2. * M_PI, 3);
  // Assemble the terms correctly to obtain the weighted densities
  for (size_t i = 0; i != grid_count; ++i) {
    r = dr * static_cast<double>(i+1);
    rr = r * r;
    rrr = rr * r;
    scalar_weighted_dens_real->at(0, i) =
        scalar_weighted_dens_four->at(0, i+1) / r;
    scalar_weighted_dens_real->at(1, i) =
        scalar_weighted_dens_four->at(1, i+1) / r;
    scalar_weighted_dens_real->at(2, i) =
        scalar_weighted_dens_four->at(2, i+1) / r;
    scalar_weighted_dens_real->at(3, i) =
        scalar_weighted_dens_four->at(3, i+1) / r;
    vector_weighted_dens_real->at(0, i) =
        vector_weighted_dens_four->at(0, i+1) / rr -
        vector_weighted_dens_four->at(1, i+1) / r;
    vector_weighted_dens_real->at(1, i) =
        vector_weighted_dens_four->at(2, i+1) / rr -
        vector_weighted_dens_four->at(3, i+1) / r;
    tensor_weighted_dens_real->at(0, i) =
        tensor_weighted_dens_four->at(0, i+1) / rrr -
        tensor_weighted_dens_four->at(1, i+1) / rr -
        tensor_weighted_dens_four->at(2, i+1) / (3. * r);
    tensor_weighted_dens_real->at(1, i) =
        2. * tensor_weighted_dens_four->at(1, i+1) / rr +
        2. * tensor_weighted_dens_four->at(2, i+1) / (3. * r) -
        2. * tensor_weighted_dens_four->at(0, i+1) / rrr;
  }
  // The convolution may have numerical problems at the outermost part.
  // One way of avoiding this is by using an outer external potential
  // forcing the density profile to vanish (i.e. sperical capacitor).
}
// _____________________________________________________________________________
void FunctionalFMTSpherical::check_weighted_densities() {
  // initialize counters
  size_t n3bigger1_counter = 0;
  size_t n3smaller0_counter = 0;
  size_t n2smaller0_counter = 0;
  // Since calc_weighted_densities() is executed before this function,
  // we can check if the weighted densities are in the legal range.
  for (size_t i = 0; i < grid_count; ++i) {
    // Check scalar weighted densities for values smaller 0
    if (scalar_weighted_dens_real->at(0, i) < 0.) {
      ++n3smaller0_counter;
      scalar_weighted_dens_real->at(0, i) = 0.;
    }
    if (scalar_weighted_dens_real->at(1, i) < 0.) {
      ++n2smaller0_counter;
      scalar_weighted_dens_real->at(1, i) = 0.;
      scalar_weighted_dens_real->at(2, i) = 0.;
      scalar_weighted_dens_real->at(3, i) = 0.;
    }
    // Check if n3 is bigger than 1
    if (scalar_weighted_dens_real->at(0, i) > 1.) {
      ++n3bigger1_counter;
      scalar_weighted_dens_real->at(0, i) = 1.;
    }
  }
  // Show warning
  if (n3smaller0_counter > 0) {
    std::cerr << "FunctionalFMTSpherical::check_weighted_densities(): \"";
    std::cerr << "Warning: Local n3 < 0.0 at " << n3smaller0_counter;
    std::cerr << " positions.\"" << std::endl;
  }
  if (n2smaller0_counter > 0) {
    std::cerr << "FunctionalFMTSpherical::check_weighted_densities(): \"";
    std::cerr << "Warning: Local n2 < 0.0 at " << n2smaller0_counter;
    std::cerr << " positions.\"" << std::endl;
  }
  if (n3bigger1_counter > 0) {
    std::cerr << "FunctionalFMTSpherical::check_weighted_densities(): \"";
    std::cerr << "Warning: Local n3 > 1.0 at " << n3bigger1_counter;
    std::cerr << " positions.\"" << std::endl;
  }
}
// _____________________________________________________________________________
void FunctionalFMTSpherical::calc_partial_derivatives() {
  // Calculate the partial deriavative at every position
  for (size_t i = 0; i != grid_count; ++i) {
    calc_local_partial_derivatives(i);
  }
}
// _____________________________________________________________________________
void FunctionalFMTSpherical::calc_local_partial_derivatives(size_t i) {
  // Auxiliary variables for the weight functions
  // Skalars
  double n3, n3n3, n3n3n3, oneMn3, oneOoneMn3, oneOn3, oneOn3n3;
  double n2, n2n2;
  double n1;
  double n0;
  // Vectors: only radial component counts
  double nvec2, nvec2nvec2;
  double nvec1;
  // Tensor: diagonal, 3 radial components
  double ntensorm2first, ntensorm2second, ntensorm2third, trace3;
  // Auxiliary variables \phi^{num}_2 and \phi^{num}_3 (B.18)-(B.21)
  double phi2, phi3;
  // Auxiliary variables \partial\phi^{num}_j/partial n_3 (B.30)-(B.33)
  double dphi2dn3, dphi3dn3;
  // Auxiliary variables \partial\Phi_j/\partial n_k (B.37)-(B.46)
  double dPhi1dn0, dPhi1dn3;
  double dPhi2dn3, dPhi2dn2, dPhi2dn1, dPhi2dnvec2, dPhi2dnvec1;
  double dPhi3dn2, dPhi3dn3, dPhi3dnvec2;
  double dPhi3dnmat2first, dPhi3dnmat2third;
  // Auxiliary constants
  double oneO24pi;
  // Calculate auxiliary weight functions ("M"=minus, "O"=over)
  n3 = scalar_weighted_dens_real->at(0, i);  // n3
  n3n3 = n3 * n3;
  n3n3n3 = n3n3 * n3;
  oneMn3 = 1. - n3;
  oneOoneMn3 = 1. / oneMn3;
  oneOn3 = 1. / n3;
  oneOn3n3 = 1. / n3n3;
  n2 = scalar_weighted_dens_real->at(1, i);  // n2
  n2n2 = n2 * n2;
  n1 = scalar_weighted_dens_real->at(2, i);  // n1
  n0 = scalar_weighted_dens_real->at(3, i);  // n0
  nvec2 = vector_weighted_dens_real->at(0, i);  // nvec2
  nvec2nvec2 = nvec2 * nvec2;
  nvec1 = vector_weighted_dens_real->at(1, i);  // nvec1
  ntensorm2first = tensor_weighted_dens_real->at(0, i);  // tensor
  ntensorm2second = ntensorm2first;
  ntensorm2third = tensor_weighted_dens_real->at(1, i);  // tensor
  trace3 = ntensorm2first * ntensorm2first * ntensorm2first +
      ntensorm2second * ntensorm2second * ntensorm2second  +
      ntensorm2third * ntensorm2third * ntensorm2third;
  // Calculate auxiliary constants
  oneO24pi = 1. / (24. * M_PI);
  // If argument of log is close to one, use Taylor series, because of 1/n3.
  if (n3 >= 1.0e-6) {
    phi2 = (5./3.) + (2./3.) * (oneMn3*oneOn3) * log(oneMn3) - (n3/3.);
    phi3 = 2. - (2./3.) *
        (oneOn3 + n3 + (oneMn3 * oneMn3 * oneOn3n3) * log(oneMn3));
    dphi2dn3 = -(2./3.) * (oneOn3n3 * log(oneMn3) + oneOn3 + .5);
    dphi3dn3 = (2./3.) * oneOn3n3 *
        (2. - n3 - n3n3 + 2. * (oneMn3*oneOn3) * log(oneMn3));
  } else {
    phi2 = 1. + (1./9.) * n3n3 + (1./18.) * n3n3n3 + (1./30.) * n3n3n3 * n3 +
        (1./45.) * n3n3n3 * n3n3;
    phi3 = 1. - (4./9.) * n3 + (1./18.) * n3n3 + (1./45.) * n3n3n3 +
        (1./90.) * n3n3n3 * n3;
    dphi2dn3 = (2./9.) * n3 + (1./6.) * n3n3 + (2./15) * n3n3n3 +
        (1./9.) * n3n3n3 * n3;
    dphi3dn3 = -(4./9.) + (1./9.) * n3 + (1./15.) * n3n3 + (2./45.) * n3n3n3;
  }
  // Calculate the partial derivatives \partial\Phi_j/\partial n_k (B.37)-(B.46)
  // Phi1 derivatives
  dPhi1dn3 = n0 * oneOoneMn3;
  dPhi1dn0 = -log(oneMn3);
  // Phi2 derivatives
  dPhi2dn3 = ((phi2*oneOoneMn3) + dphi2dn3) * (n2 * n1 - nvec1 * nvec2) *
      oneOoneMn3;
  dPhi2dn2 = phi2 * n1 * oneOoneMn3;
  dPhi2dn1 = phi2 * n2 * oneOoneMn3;
  dPhi2dnvec2 = -phi2 * nvec1 * oneOoneMn3;  // times direction vector
  dPhi2dnvec1 = -phi2 * nvec2 * oneOoneMn3;  // times direction vector
  // Phi3 derivatives
  dPhi3dn2 = phi3 * (3. * n2n2 - 3. * nvec2nvec2) * oneO24pi * oneOoneMn3 *
      oneOoneMn3;
  dPhi3dn3 = ((2. * phi3 * oneOoneMn3) + dphi3dn3) *
      (n2n2 * n2 - 3. * n2 * nvec2nvec2 +
      4.5 * (nvec2nvec2 * ntensorm2third - trace3)) * oneO24pi * oneOoneMn3 *
      oneOoneMn3;
  dPhi3dnvec2 = phi3 * (9. * ntensorm2third * nvec2 - 6. * n2 * nvec2) *
      oneO24pi * oneOoneMn3 * oneOoneMn3;  // times direction vector
  dPhi3dnmat2first = 4.5 * phi3 * (-3. * ntensorm2first * ntensorm2first) *
      oneO24pi * oneOoneMn3 * oneOoneMn3;  // times direction tensor
  dPhi3dnmat2third = 4.5 * phi3 *
      (nvec2nvec2 - 3. * ntensorm2third * ntensorm2third) *
      oneO24pi * oneOoneMn3 * oneOoneMn3;  // times direction tensor
  // Sum partial derivatives \sum_{i=1}^3\partial\Phi_i/\partial n_k
  // Scalar n3
  scalar_partial_derivative_real->at(0, i) = dPhi1dn3 + dPhi2dn3 + dPhi3dn3;
  // Scalar n2
  scalar_partial_derivative_real->at(1, i) = dPhi2dn2 + dPhi3dn2;
  // Scalar n1
  scalar_partial_derivative_real->at(2, i) = dPhi2dn1;
  // Scalar n0
  scalar_partial_derivative_real->at(3, i) = dPhi1dn0;
  // Vector n2
  vector_partial_derivative_real->at(0, i) = dPhi2dnvec2 + dPhi3dnvec2;
  // Vector n1
  vector_partial_derivative_real->at(1, i) = dPhi2dnvec1;
  // Tensor nm2 (first and second diagonal element are the same)
  tensor_partial_derivative_real->at(0, i) = dPhi3dnmat2first;
  tensor_partial_derivative_real->at(1, i) = dPhi3dnmat2third;
}
// _____________________________________________________________________________
void FunctionalFMTSpherical::calc_weighted_partial_derivatives(
    DataField<double>* functional_derivative) {
  // TODO(Moritz): speed anhancement by not using "at" and by saving transforms
  double r{0.};
  double kr{0.}, krkr{0.};
  std::vector<fftw_plan> forward_plans;
  std::vector<fftw_plan> backward_plans;
  // Specify the plans
  for (size_t i = 0; i != 4; ++i) {  // TODO(Moritz): make number static
    forward_plans.push_back(
        fftw_plan_r2r_1d(grid_count,
            scalar_derivative_terms_four->array(i) + 1,
            scalar_derivative_terms_four->array(i) + 1,
            FFTW_RODFT00, flags_destroy));
  }
  forward_plans.push_back(
      fftw_plan_r2r_1d(grid_count,
          vector_derivative_terms_four->array(0) + 1,
          vector_derivative_terms_four->array(0) + 1,
          FFTW_RODFT00, flags_destroy));
  forward_plans.push_back(
      fftw_plan_r2r_1d(grid_count + 1,
          vector_derivative_terms_four->array(1),
          vector_derivative_terms_four->array(1),
          FFTW_REDFT00, flags_destroy));
  forward_plans.push_back(
      fftw_plan_r2r_1d(grid_count,
          vector_derivative_terms_four->array(2) + 1,
          vector_derivative_terms_four->array(2) + 1,
          FFTW_RODFT00, flags_destroy));
  forward_plans.push_back(
      fftw_plan_r2r_1d(grid_count + 1,
          vector_derivative_terms_four->array(3),
          vector_derivative_terms_four->array(3),
          FFTW_REDFT00, flags_destroy));
  forward_plans.push_back(
      fftw_plan_r2r_1d(grid_count,
          tensor_derivative_terms_four->array(0) + 1,
          tensor_derivative_terms_four->array(0) + 1,
          FFTW_RODFT00, flags_destroy));
  forward_plans.push_back(
      fftw_plan_r2r_1d(grid_count,
          tensor_derivative_terms_four->array(1) + 1,
          tensor_derivative_terms_four->array(1) + 1,
          FFTW_RODFT00, flags_destroy));
  forward_plans.push_back(
      fftw_plan_r2r_1d(grid_count + 1,
          tensor_derivative_terms_four->array(2),
          tensor_derivative_terms_four->array(2),
          FFTW_REDFT00, flags_destroy));
  forward_plans.push_back(
      fftw_plan_r2r_1d(grid_count,
          tensor_derivative_terms_four->array(3) + 1,
          tensor_derivative_terms_four->array(3) + 1,
          FFTW_RODFT00, flags_destroy));
  for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
    backward_plans.push_back(
        fftw_plan_r2r_1d(grid_count,
            functional_derivative->array(*it),
            functional_derivative->array(*it),
            FFTW_RODFT00, flags_destroy));
  }
  // Transform the partial derivatives into Fourier space.
  // The vectorial and tensorial derivatives need again a special transform.
  // The dummy arrays dfdnrs_terms and dfdns_terms are used.
  for (size_t i = 0; i != 4; ++i) {
    scalar_derivative_terms_four->at(i, 0) = 0.;
    vector_derivative_terms_four->at(i, 0) = 0.;
    tensor_derivative_terms_four->at(i, 0) = 0.;
  }
  for (size_t i = 1; i < grid_count + 1; ++i) {
    r = static_cast<double>(i) * dr;
    // Partial derivatives of scalar weighted densities
    scalar_derivative_terms_four->at(0, i) =
        scalar_partial_derivative_real->at(0, i-1) * r;  // sin*r of n3
    scalar_derivative_terms_four->at(1, i) =
        scalar_partial_derivative_real->at(1, i-1) * r;  // sin*r of n2
    scalar_derivative_terms_four->at(2, i) =
        scalar_partial_derivative_real->at(2, i-1) * r;  // sin*r of n1
    scalar_derivative_terms_four->at(3, i) =
        scalar_partial_derivative_real->at(3, i-1) * r;  // sin*r of n0
    // Partial derivatives of vector weighted densities
    vector_derivative_terms_four->at(0, i) =
        vector_partial_derivative_real->at(0, i-1);  // sin of nvec2
    vector_derivative_terms_four->at(1, i) =
        vector_partial_derivative_real->at(0, i-1) * r;  // cos*r of nvec2
    vector_derivative_terms_four->at(2, i) =
        vector_partial_derivative_real->at(1, i-1);  // sin of nvec1
    vector_derivative_terms_four->at(3, i) =
        vector_partial_derivative_real->at(1, i-1) * r;  // cos*r of nvec1
    // Partial derivatives of tensor weighted densities
    tensor_derivative_terms_four->at(0, i) =
        tensor_partial_derivative_real->at(0, i-1) * r;
    tensor_derivative_terms_four->at(1, i) =
        (tensor_partial_derivative_real->at(1, i-1) -
        tensor_partial_derivative_real->at(0, i-1)) / r;
    tensor_derivative_terms_four->at(2, i) =
        (tensor_partial_derivative_real->at(1, i-1) -
        tensor_partial_derivative_real->at(0, i-1));
    tensor_derivative_terms_four->at(3, i) =
        (tensor_partial_derivative_real->at(1, i-1) -
        tensor_partial_derivative_real->at(0, i-1)) * r;
  }
  for (auto it = forward_plans.begin(); it != forward_plans.end(); ++it) {
    fftw_execute(*it);
  }
  // Just to be safe
  for (size_t i = 0; i != 4; ++i) {
    scalar_derivative_terms_four->at(i, 0) = 0.;
    vector_derivative_terms_four->at(i, 0) = 0.;
    tensor_derivative_terms_four->at(i, 0) = 0.;
  }
  // Normalize
  for (size_t i = 1; i != grid_count+1; ++i) {
    kr = dkr * static_cast<double>(i);
    krkr = kr * kr;
    for (size_t j = 0; j != 4; ++j) {
      scalar_derivative_terms_four->at(j, i) *= 4. * M_PI * norm_sin / kr;
    }
    vector_derivative_terms_four->at(0, i) *= 4. * M_PI * norm_sin / krkr;
    vector_derivative_terms_four->at(1, i) *= 4. * M_PI * norm_cos / kr;
    vector_derivative_terms_four->at(2, i) *= 4. * M_PI * norm_sin / krkr;
    vector_derivative_terms_four->at(3, i) *= 4. * M_PI * norm_cos / kr;
    tensor_derivative_terms_four->at(0, i) *= 4. * M_PI * norm_sin / kr;
    tensor_derivative_terms_four->at(1, i) *= 4. * M_PI * norm_sin /(krkr * kr);
    tensor_derivative_terms_four->at(2, i) *= 4. * M_PI * norm_cos / krkr;
    tensor_derivative_terms_four->at(3, i) *= 4. * M_PI * norm_sin / kr;
  }
  // Merge terms
  for (size_t i = 0; i < grid_count; ++i) {
    scalar_derivative_four->at(0, i) = scalar_derivative_terms_four->at(0, i+1);
    scalar_derivative_four->at(1, i) = scalar_derivative_terms_four->at(1, i+1);
    scalar_derivative_four->at(2, i) = scalar_derivative_terms_four->at(2, i+1);
    scalar_derivative_four->at(3, i) = scalar_derivative_terms_four->at(3, i+1);
    vector_derivative_four->at(0, i) =
        vector_derivative_terms_four->at(1, i+1) -
        vector_derivative_terms_four->at(0, i+1);  // nvec2
    vector_derivative_four->at(1, i) =
        vector_derivative_terms_four->at(3, i+1) -
        vector_derivative_terms_four->at(2, i+1);  // nvec1
    tensor_derivative_four->at(0, i) =
        tensor_derivative_terms_four->at(0, i+1) +
        tensor_derivative_terms_four->at(1, i+1) -
        tensor_derivative_terms_four->at(2, i+1);  // 1st and 2nd el. of tensor
    tensor_derivative_four->at(1, i) =
        tensor_derivative_terms_four->at(0, i+1) -
        2. * tensor_derivative_terms_four->at(1, i+1) +
        2. * tensor_derivative_terms_four->at(2, i+1) +
        tensor_derivative_terms_four->at(3, i+1);  // 3rd el. of tensor
  }
  functional_derivative->zeros();
  size_t spec_i{0};
  for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
    spec_i = it - affected_species.begin();
    for (size_t j = 0; j != grid_count; ++j) {
      kr = dkr * static_cast<double>(j + 1);
      // Scalar
      functional_derivative->at(*it, j) += scalar_derivative_four->at(0, j) *
          kr * weights_four.at(spec_i).at(0, j+1);
          // hardcoded sign (sym.:positive)
      functional_derivative->at(*it, j) += scalar_derivative_four->at(1, j) *
          kr * weights_four.at(spec_i).at(1, j+1);
          // hardcoded sign (sym.:positive)
      functional_derivative->at(*it, j) += scalar_derivative_four->at(2, j) *
          kr * weights_four.at(spec_i).at(1, j+1) /
          (2. * M_PI * diameters.at(spec_i));
          // hardcoded sign (sym.:positive)
      functional_derivative->at(*it, j) += scalar_derivative_four->at(3, j) *
          kr * weights_four.at(spec_i).at(1, j+1) /
          (M_PI * diameters.at(spec_i) * diameters.at(spec_i));
          // hardcoded sign (sym.:positive)
      // Vector
      functional_derivative->at(*it, j) += -vector_derivative_four->at(0, j) *
          kr * kr * weights_four.at(spec_i).at(0, j+1);
          // hardcoded sign (asym.:negative)
      functional_derivative->at(*it, j) += -vector_derivative_four->at(1, j) *
          kr * kr * weights_four.at(spec_i).at(0, j+1) /
          (2. * M_PI * diameters.at(spec_i));
          // hardcoded sign (asym.:negative)
      // Tensor
      functional_derivative->at(*it, j) += (2. / 3.) *
          (tensor_derivative_four->at(1, j) -
          tensor_derivative_four->at(0, j)) * kr *
          weights_four.at(spec_i).at(2, j+1);  // hardcoded sign (sym.:positive)
    }
  }
  // Back transform
  for (auto it = backward_plans.begin(); it != backward_plans.end(); ++it) {
    fftw_execute(*it);
  }
  for (size_t i = 0; i != grid_count; ++i) {
    r = static_cast<double>(i + 1) * dr;
    for (auto it = affected_species.begin(); it != affected_species.end();
        ++it) {
      functional_derivative->at(*it, i) *= 4. * M_PI / r;
    }
  }
  *(functional_derivative) *= norm_sin;
  // 3D FT Backtransform normalization
  *(functional_derivative) *= 1. / pow(2. * M_PI, 3);
}
// _____________________________________________________________________________
double FunctionalFMTSpherical::calc_local_energy_density(size_t position) {
  // Define a lot of dummy variables
  double Phi1, Phi2, Phi3;  // energy density terms
  double phi2, phi3;  // factors in the energy density terms
  double phi2_num, phi3_num;  // factors in the energy density terms
  double n0, n1, n2, n3;  // scalar weighted densities
  double n3n3, n3n3n3, oneMn3, logOneMn3, oneMn3squared;  // auxiliary variables
  double n2n2, n2n2n2;  // auxiliary variables
  double trace2, trace3;  // auxiliary variables
  double nvec1, nvec2;  // vectorial weighted densities (only z-comp. is non-0)
  double ntensorm2first;  // tensorial weighted densities (tensor is diagonal)
  double ntensorm2second;
  double ntensorm2third;
  // Assign values to the weighted densities
  n3 = scalar_weighted_dens_real->at(0, position);
  n2 = scalar_weighted_dens_real->at(1, position);
  n1 = scalar_weighted_dens_real->at(2, position);
  n0 = scalar_weighted_dens_real->at(3, position);
  nvec2 = vector_weighted_dens_real->at(0, position);
  nvec1 = vector_weighted_dens_real->at(1, position);
  ntensorm2first = tensor_weighted_dens_real->at(0, position);
  ntensorm2second = tensor_weighted_dens_real->at(0, position);
  ntensorm2third = tensor_weighted_dens_real->at(1, position);
  // Calculate auxiliary variables
  n3n3 = n3 * n3;
  n3n3n3 = n3 * n3 * n3;
  oneMn3 = 1. - n3;
  logOneMn3 = log(oneMn3);
  oneMn3squared = oneMn3 * oneMn3;
  n2n2 = n2 * n2;
  n2n2n2 = n2 * n2n2;
  trace2 = ntensorm2first * ntensorm2first + ntensorm2second * ntensorm2second +
      ntensorm2third * ntensorm2third;
  trace3 = ntensorm2first * ntensorm2first * ntensorm2first +
      ntensorm2second * ntensorm2second * ntensorm2second +
      ntensorm2third * ntensorm2third * ntensorm2third;
  // Calculate factors in the energy density terms
  if (n3 < 1e-5) {  // avoiding logarithm of very small numbers
    phi2 = 1. + .5 * n3 + .3 * n3n3 + .2 * n3n3n3;  // +O(n^4)
    phi3 = 1. - .125 * n3 - .05 * n3n3 - .025 * n3n3n3;  // +O(n^4)
  } else {
    phi2 = (6. * n3 - 3. * n3n3 + 6. * oneMn3 * logOneMn3) / n3n3n3;
    phi3 = (6. * n3 - 9. * n3n3 + 6. * n3n3n3 + 6. * oneMn3squared * logOneMn3)/
        (4. * n3n3n3);
  }
  phi2_num = 1. + n3n3 * phi2 / 9.;
  phi3_num = 1. - 4. * n3 * phi3 / 9.;
  // Calculate free energy density terms
  Phi1 = -n0 * logOneMn3;
  Phi2 = phi2_num *  (n1 * n2 - nvec1 * nvec2) / oneMn3;
  Phi3 = phi3_num * (n2n2n2 - 3. * n2 * nvec2 * nvec2 + 4.5 * (
      n2n2 * ntensorm2third - n2 * nvec2 * nvec2 - trace3 + n2 * trace2)) /
      (24. * M_PI * oneMn3squared);
  return Phi1 + Phi2 + Phi3;
}
// _____________________________________________________________________________
double FunctionalFMTSpherical::radial_integration(
    double* data, int n, double delta) {
  // Integrate with closed Newton-Cotes formula: Num. Rep. 3rd ed. eq. 4.1.14.
  // r = 0 has no contribution.
  double integral = 0.;
  double r = 0.;
  if (grid_count < 6) {
    std::cerr << "FunctionalFMTSpherical::radial_integration(): ";
    std::cerr << "\"Error: Integration needs more grid points.\"";
    std::cerr << std::endl;
    return 0.;
  }
  for (size_t i = 0; i != grid_count; ++i) {
    r = static_cast<double>(i + 1) * dr;
    if (i == grid_count - 1)
      integral += data[i] * r * r * (3. / 8.);
    else if (i == 0 || i == grid_count - 2)
      integral += data[i] * r * r * (7. / 6.);
    else if (i == 1 || i == grid_count - 3)
      integral += data[i] * r * r * (23. / 24.);
    else
      integral += data[i] * r * r;
  }
  // Spherical coordinates --> 4\pi
  integral *= delta * 4. * M_PI;
  return integral;
}
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
