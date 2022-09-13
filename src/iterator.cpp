// SPDX-FileCopyrightText: 2021 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "iterator.hpp"  // NOLINT
#include <cmath>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <limits>
#include "data_frame.hpp"  // NOLINT
#include "convergence_criterion.hpp"  // NOLINT
#include "convergence_criterion_max_dev.hpp"  // NOLINT
#include "convergence_criterion_steps.hpp"  // NOLINT
#include "convergence_criterion_nan.hpp"  // NOLINT
// _____________________________________________________________________________
Iterator::Iterator(
    std::vector<DataFrame<1, double>>* density_profiles,
    const std::vector<DataFrame<1, double>>& exp_external_potentials,
    const std::vector<Properties>& species_properties)
  : density_profiles(density_profiles),
    exp_external_potentials(&exp_external_potentials),
    species_properties(&species_properties) {
  size_t grid_count = density_profiles->at(0).size();
  size_t species_count = density_profiles->size();
  proposed_densities = std::vector(species_count,
      DataFrame<1, double>(grid_count));
  add_convergence_criterion<ConvergenceCriterionMaxDev>(1.0e-6);
  add_convergence_criterion<ConvergenceCriterionSteps>(1e4);
  add_convergence_criterion<ConvergenceCriterionNan>(0);
}
// _____________________________________________________________________________
Iterator::~Iterator() {
  clear_convergence_criteria();
}
// _____________________________________________________________________________
void Iterator::add_excess_functional(Functional* functional) {
  size_t species_count{species_properties->size()};
  size_t index{0};
  excess_functionals.push_back(functional);
  // For every functional and every species a DataFrame is created
  bulk_derivatives.push_back(std::vector<double>(0));
  functional_derivatives.push_back(std::vector<DataFrame<1, double>>(0));
  index = functional_derivatives.size() - 1;
  for (size_t i = 0; i < species_count; ++i) {
    bulk_derivatives.at(index).push_back(0.);
    functional_derivatives.at(index).push_back(DataFrame<1, double>(
        density_profiles->at(0).size()));
  }
}
// _____________________________________________________________________________
void Iterator::remove_excess_functional(size_t index) {
  excess_functionals.erase(excess_functionals.begin() + index);
  functional_derivatives.erase(functional_derivatives.begin() + index);
  bulk_derivatives.erase(bulk_derivatives.begin() + index);
}
// _____________________________________________________________________________
void Iterator::clear_functionals() {
  excess_functionals.clear();
  functional_derivatives.clear();
  bulk_derivatives.clear();
}
// _____________________________________________________________________________
void Iterator::run_picard(double mixing) {
  double bulk_density{0.};
  size_t steps{0};
  size_t species_count{species_properties->size()};
  // Calculate bulk derivatives
  for (size_t i = 0; i < excess_functionals.size(); ++i) {
    excess_functionals.at(i)->calc_bulk_derivative(&bulk_derivatives.at(i));
  }
  // Iterator loop
  while (!check_convergence_criteria()) {
    ++steps;
    // Calculate the functional derivatives.
    for (size_t i = 0; i < excess_functionals.size(); ++i) {
      excess_functionals.at(i)->calc_derivative(&functional_derivatives.at(i));
    }
    // Picard iteration scheme
    for (size_t i = 0; i < species_count; ++i) {
      species_properties->at(i).get_property("bulk density", &bulk_density);
      // Calculate the right hand side of the update formula of cDFT.
      proposed_densities.at(i) = bulk_density * exp_external_potentials->at(i);
      for (size_t j = 0; j < functional_derivatives.size(); ++j) {
        proposed_densities.at(i) *= exp(bulk_derivatives.at(j).at(i)) *
            exp(-1. * functional_derivatives.at(j).at(i));
      }
      // Mix the new solution with the old one
      density_profiles->at(i) = density_profiles->at(i) * (1. - mixing) +
          mixing * proposed_densities.at(i);
    }
  }
}
// _____________________________________________________________________________
void Iterator::run_anderson(double mixing, size_t memory) {
  double bulk_density{0.};
  size_t steps{0};
  size_t species_count{species_properties->size()};
  size_t grid_count{proposed_densities.at(0).size()};
  // A lot of memory is needed for all the past profiles
  std::vector<DataFrame<1, double>> deviation_profiles(
      species_count, DataFrame<1, double>(grid_count));
  std::vector<std::vector<DataFrame<1, double>>> past_density_profiles(0);
  std::vector<std::vector<DataFrame<1, double>>> past_proposed_densities(0);
  std::vector<std::vector<DataFrame<1, double>>> past_density_deviations(0);
  std::vector<std::vector<double>> scalar_products(0);
  std::vector<double> alphas(0);
  // Calculate bulk derivatives
  for (size_t i = 0; i < excess_functionals.size(); ++i) {
    excess_functionals.at(i)->calc_bulk_derivative(&bulk_derivatives.at(i));
  }
  // Iterator loop
  while (!check_convergence_criteria()) {
    ++steps;
    // Calculate the functional derivatives.
    for (size_t i = 0; i < excess_functionals.size(); ++i) {
      excess_functionals.at(i)->calc_derivative(&functional_derivatives.at(i));
    }
    // Calculate new density proposal via the Picard scheme
    for (size_t i = 0; i < species_count; ++i) {
      species_properties->at(i).get_property("bulk density", &bulk_density);
      // Calculate the right hand side of the update formula of cDFT.
      proposed_densities.at(i) = bulk_density * exp_external_potentials->at(i);
      for (size_t j = 0; j < functional_derivatives.size(); ++j) {
        proposed_densities.at(i) *= exp(bulk_derivatives.at(j).at(i)) *
            exp(-1. * functional_derivatives.at(j).at(i));
        deviation_profiles.at(i) =
            proposed_densities.at(i) - density_profiles->at(i);
      }
    }
    // Append density profile, proposed density profile and deviation profile
    // to history
    past_density_profiles.push_back(*(density_profiles));
    past_proposed_densities.push_back(proposed_densities);
    past_density_deviations.push_back(deviation_profiles);
    if (past_density_profiles.size() > memory) {
      past_density_profiles.erase(past_density_profiles.begin());
      past_proposed_densities.erase(past_proposed_densities.begin());
      past_density_deviations.erase(past_density_deviations.begin());
    }
    // Calculate all possible scalar products between the density deviations
    scalar_products = calc_scalar_products(past_density_deviations);
    // Find the linear combination (alphas) that minimizes the norm
    alphas = shortest_linear_combination(scalar_products);
    // Update formula
    for (size_t i = 0; i < species_count; ++i) {
      density_profiles->at(i).zero();
      for (size_t j = 0; j < past_density_profiles.size(); ++j) {
        density_profiles->at(i) += (1. - mixing) * alphas.at(j) *
            past_density_profiles.at(j).at(i);
        density_profiles->at(i) += mixing * alphas.at(j) *
            past_proposed_densities.at(j).at(i);
      }
    }
  }
}
// _____________________________________________________________________________
double Iterator::calculate_excess_free_energy() {
  double energy{0.};
  for (auto functional : excess_functionals) {
    energy += functional->calc_energy();
  }
  return energy;
}
// _____________________________________________________________________________
double Iterator::calculate_gc_energy() {
  double energy{0.};
  // Calculate excess free energy
  energy += calculate_excess_free_energy();
  // Calculate ideal free energy
  // TODO(Moritz):
  // Calculate external and chemical potential term
  // TODO(Moritz):
  return energy;
}
// _____________________________________________________________________________
void Iterator::clear_convergence_criteria() {
  for (auto criterion : convergence_criteria) {
    delete(criterion);
  }
  convergence_criteria.clear();
}
// _____________________________________________________________________________
bool Iterator::check_convergence_criteria() {
  bool converged{false};
  double progress{0.};
  std::cout << "Iterator::check_convergence_criteria(): \"";
  for (auto criterion : convergence_criteria) {
    converged = (converged | criterion->check(&progress));
    std::cout << criterion->name() << ": " << progress << "; ";
  }
  std::cout << "\"" << std::endl;
  return converged;
}
// _____________________________________________________________________________
std::vector<std::vector<double>> Iterator::calc_scalar_products(
    std::vector<std::vector<DataFrame<1, double>>>& profile_history) {
  size_t hist_size{profile_history.size()};
  size_t species_count{profile_history.at(0).size()};
  size_t grid_count{profile_history.at(0).at(0).size()};
  std::vector<std::vector<double>> scalar_products(
      hist_size, std::vector<double>(hist_size));
  double sum {0.};
  for (size_t i = 0; i < hist_size; ++i) {
    for (size_t j = i; j < hist_size; ++j) {
      sum = 0.;
      for (size_t k = 0; k < species_count; ++k) {
        for (size_t l = 0; l < grid_count; ++l) {
          sum += profile_history.at(i).at(k).at(l) * profile_history.at(j).at(k).at(l);
        }
      }
      scalar_products.at(i).at(j) = scalar_products.at(j).at(i) = sum;
    }
  }
  return scalar_products;
}
// _____________________________________________________________________________
std::vector<double> Iterator::shortest_linear_combination(
    std::vector<std::vector<double>>& scalar_products) {
  size_t vector_count{scalar_products.size()};
  size_t param_count{vector_count + 1};  // Lagrange multiplier gives +1
  size_t iter{0};
  int status{GSL_CONTINUE};
  std::vector<double> alphas(vector_count);
  // Initial guess
  gsl_vector* x = gsl_vector_calloc(param_count);
  for (size_t i = 0; i < vector_count; ++i) {
    gsl_vector_set(x, i, 1. / static_cast<double>(vector_count));
  }
  // Find the root
  const gsl_multiroot_fsolver_type* T = gsl_multiroot_fsolver_hybrids;
  //const gsl_multiroot_fsolver_type* T = gsl_multiroot_fsolver_dnewton;
  gsl_multiroot_fsolver* s = gsl_multiroot_fsolver_alloc (T, param_count);
  gsl_multiroot_function f = {&anderson_f, param_count, &scalar_products};
  gsl_multiroot_fsolver_set (s, &f, x);
  while ((status == GSL_CONTINUE) && (iter++ < 1000)) {
    status = gsl_multiroot_fsolver_iterate(s);
    status = gsl_multiroot_test_residual (s->f, 1e-10);
  }
  for (size_t i = 0; i < vector_count; ++i) {
    alphas.at(i) = gsl_vector_get(s->x, i);
  }
  return alphas;
}
// _____________________________________________________________________________
int anderson_f (const gsl_vector* x, void* params, gsl_vector* f) {
  std::vector<std::vector<double>> scalar_products =
      *((std::vector<std::vector<double>>*) params);
  size_t vector_count{scalar_products.size()};
  double sum{0.};
  double lambda{gsl_vector_get(x, vector_count)};
  for (size_t i = 0; i < vector_count; ++i) {
    sum = 0.;
    for (size_t j = 0; j < vector_count; ++j) {
      sum += gsl_vector_get(x, j) * scalar_products.at(j).at(i);
    }
    gsl_vector_set(f, i, 2. * sum - lambda);
  }
  sum = 0.;
  for (size_t i = 0; i < vector_count; ++i) {
    sum += gsl_vector_get(x, i);
  }
  gsl_vector_set(f, vector_count, sum - 1.);
  return GSL_SUCCESS;
}
// _____________________________________________________________________________
// _____________________________________________________________________________
