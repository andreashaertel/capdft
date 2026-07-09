// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file boundaries.cpp
 *  \brief Source file for the Boundaries class.
 *
 *  The file contains the definitions of the Boundaries class.
 */
#include <list>
#include "boundaries.hpp"
#include "boundary_surface.hpp"
#include "system.hpp"
#include "discretization_handler.hpp"

// Template class instantiation
template class Boundaries<Cartesian>;

// _____________________________________________________________________________
template <Geometry geometry>
Boundaries<geometry>::Boundaries() {
}
// _____________________________________________________________________________
template <Geometry geometry>
void Boundaries<geometry>::add_surface(BoundarySurface<geometry>* surface) {
  this->push_back(surface);
}  
// _____________________________________________________________________________
template <Geometry geometry>
void Boundaries<geometry>::set_all_boundary_values(double value) {
  for (auto i = this->begin(); i != this->end(); i++) {
    (*i)->set_boundary_value(value);
  }
}
// _____________________________________________________________________________
template <Geometry geometry>
void Boundaries<geometry>::exp_external_potential_hs(System<geometry>& system,
		std::vector<DataFrame<dimension(geometry), double>>* result,
		double resolution) {
  std::vector<size_t> grid_counts = result->at(0).size_dim();
  // if there are no boundary objects, set all to 1.
  if (this->size() == 0) {
    DataFrame<dimension(geometry), double> dummy(grid_counts);
    dummy.set_all_elements_to(1.);
    for (size_t i = 0; i < result->size(); i++) {
      result->at(i) = dummy;
    }
    return;
  } else {
    auto it = this->begin();
    // initialize with first boundary's potential
    if (resolution > 0.) {
      (*it)->exp_external_potential_hs(system, result, resolution);
    } else {
      (*it)->exp_external_potential_hs(system, result);
    }
    // multiply with all remaining boundaries' potentials
    std::vector<DataFrame<dimension(geometry), double>> dummy(result->size(), 
		    DataFrame<dimension(geometry), double>(grid_counts));
    it++;
    std::cout << "HS of first surface\n";
    while (it != this->end()) {
      if (resolution > 0.) {
        (*it)->exp_external_potential_hs(system, &dummy, resolution);
      } else {
        (*it)->exp_external_potential_hs(system, &dummy);
      }
      for (size_t j = 0; j < result->size(); j++) {
        result->at(j) *= dummy.at(j);
      }
      it++;
      std::cout << "HS of another surface\n";
    }
    return;
  }
}
// _____________________________________________________________________________
