// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "functional_fmt_spherical.hpp"
#include <fftw3.h>
#include <cmath>
#include <iostream>
// _____________________________________________________________________________
FunctionalFMTSpherical::FunctionalFMTSpherical() {
  //
}
// _____________________________________________________________________________
FunctionalFMTSpherical::FunctionalFMTSpherical(System* system) {
  // Clear all std::vectors
  diameters.clear();
  bulk_densities.clear();
  affected_species.clear();
  // Get system properties
  extract_system_properties(system);
  // Get species properties; excludes all species without diameter property
  extract_species_properties(system);
  // Initialize the density profile and density profile times radial position
  density_profile_pointer = system->get_density_profile_pointer();
  density_profile_times_r = new DataField<double>(species_count, grid_count+1);
  update_density_times_r();
  // Initialize weighted densities
  scalar_weighted_dens_real = new DataField<double>(4, grid_count);
  vector_weighted_dens_real = new DataField<double>(2, grid_count);
  tensor_weighted_dens_real = new DataField<double>(2, grid_count);
  scalar_weighted_dens_four = new DataField<double>(4, grid_count+1);
  vector_weighted_dens_four = new DataField<double>(4, grid_count+1);
  tensor_weighted_dens_four = new DataField<double>(3, grid_count+1);
  density_profile_four = new DataField<double>(species_count, grid_count+1);
  // Initialize weights
  for (size_t i = 0; i != species_count; ++i) {
    weights_four.push_back(DataField<double>(3, grid_count+1));
  }
  initialize_weights();
  // TODO(Moritz): Initialize the following objects
  // TODO(Moritz): functional_derivative
  // TODO(Moritz):
}
// _____________________________________________________________________________
FunctionalFMTSpherical::~FunctionalFMTSpherical() {
  density_profile_times_r->~DataField();
  scalar_weighted_dens_real->~DataField();
  vector_weighted_dens_real->~DataField();
  tensor_weighted_dens_real->~DataField();
  density_profile_four->~DataField();
}
// _____________________________________________________________________________
void FunctionalFMTSpherical::extract_system_properties(System* sys) {
  Properties properties(sys->get_system_properties());
  properties.get_property("length", &length);
  properties.get_property("grid count", &grid_count);
  // Calculate bin sizes
  dr = length / static_cast<double>(grid_count);
  dkr = (2. * M_PI) / (2. * static_cast<double>(grid_count+1) * dr);
}
// _____________________________________________________________________________
void FunctionalFMTSpherical::extract_species_properties(System* sys) {
  std::vector<Properties> properties = sys->get_species_properties();
  double diameter{0.};
  double bulk_density{0.};
  for (auto it = properties.begin(); it != properties.end(); ++it) {
    if (it->get_property("diameter", &diameter)) {  // only species with diam.
      if (!it->get_property("bulk density", &bulk_density)) {
        exit(1);  // TODO(Moritz): throw error
      }
      diameters.push_back(diameter);
      bulk_densities.push_back(bulk_density);
      affected_species.push_back(it - properties.begin());
    }
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
void FunctionalFMTSpherical::calc_derivative() {
  //
}
// _____________________________________________________________________________
void FunctionalFMTSpherical::calc_bulk_derivative() {
  //
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
  // TODO(Moritz): speed anhancement by not using "at" and by saving transforms
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
      density_profile_four->at(j, i) *= dr / 2.;
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
          kr * weights_four.at(j).at(0, i);
      // Scalar weighted density n_2
      scalar_weighted_dens_four->at(1, i) += density_profile_four->at(j, i) *
          kr * weights_four.at(j).at(1, i);
      // Scalar weighted density n_1
      scalar_weighted_dens_four->at(2, i) += density_profile_four->at(j, i) *
          kr * weights_four.at(j).at(1, i) / (2. * M_PI * diameters.at(j));
      // Scalar weighted density n_0
      scalar_weighted_dens_four->at(3, i) += density_profile_four->at(j, i) *
          kr * weights_four.at(j).at(1, i) / (M_PI * pow(diameters.at(j), 2));
      // Vectorial weighted density n_2 z-component in rotated system
      vector_weighted_dens_four->at(0, i) += density_profile_four->at(j, i) *
          kr * weights_four.at(j).at(0, i);  // sine transform
      vector_weighted_dens_four->at(1, i) += density_profile_four->at(j, i) *
          krkr * weights_four.at(j).at(0, i);  // cosine transform
      // Vectorial weighted density n_1 z-component in rotated system
      vector_weighted_dens_four->at(2, i) += density_profile_four->at(j, i) *
          kr * weights_four.at(j).at(0, i) /
          (2. * M_PI * diameters.at(j));  // sine transform
      vector_weighted_dens_four->at(3, i) += density_profile_four->at(j, i) *
          krkr * weights_four.at(j).at(0, i) /
          (2. * M_PI * diameters.at(j));  // cosine transform
      // Tensorial weighted density n_m2: the diagonal elements only contain
      // these three terms; all off-diagonal elements are zero
      tensor_weighted_dens_four->at(0, i) += density_profile_four->at(j, i) *
          weights_four.at(j).at(2, i) / kr;  // sin/kr
      tensor_weighted_dens_four->at(1, i) += density_profile_four->at(j, i) *
          weights_four.at(j).at(2, i);  // cos
      tensor_weighted_dens_four->at(2, i) += density_profile_four->at(j, i) *
          weights_four.at(j).at(2, i) * kr;  // sin*kr
    }
  }
  // Sine and cosine transforms
  for (auto it = backward_plans.begin(); it != backward_plans.end(); ++it) {
    fftw_execute(*it);
  }
  // Normalize
  *(scalar_weighted_dens_four) *= 4. * M_PI * dkr / 2.;
  *(vector_weighted_dens_four) *= 4. * M_PI * dkr / 2.;
  *(tensor_weighted_dens_four) *= 4. * M_PI * dkr / 2.;
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
  // TODO(Moritz):
}
// _____________________________________________________________________________
void FunctionalFMTSpherical::calc_partial_derivatives() {
  //
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
