//#include "boundary_structure.hpp"
//#include "poisson_solver.hpp"
//#include <vector>
//#include "properties.hpp"
//#include "data_frame.hpp"
//
//// TODO: this implementation is only for dim = 3
//
//// ____________________________________________________________________________
//template <size_t dim>
//void exp_external_potential_hs<dim>(BoundaryStructure<dim>& boundary_structure,
//		std::vector<Properties>& species_properties,
//		std::vector<size_t>& affected_species,
//		std::vector<DataFrame<dim, double>>* potential) {
//}
//template <>
//void exp_external_potential_hs<3>(BoundaryStructure<3>& boundary_structure,
//		std::vector<Properties>& species_properties,
//		std::vector<size_t>& affected_species,
//		std::vector<DataFrame<3, double>>* potential) {
//  // Calculate minimal distance to boundary of each point in the grid
//  std::vector<size_t> grid_counts = potential->size_dim();
//  DataFrame<dim, double> minimal_distances =
//	  boundary_structure.minimal_distances(grid_counts);
//  // Iterate over all points in the discretized grid, and compare the boundary
//  // distance to the particle radius.
//  double diameter{0.};
//  for (size_t s : affected_species) {
//    species_properties.at(s).get_property("diameter", &diameter);
//    for (size_t i = 0; i != grid_counts.at(0); ++i) {
//      for (size_t j = 0; j != grid_counts.at(1); ++j) {
//        for (size_t k = 0; k != grid_counts.at(2); ++k) {
//          if (minimal_distances.at(i, j, k) < diameter * 0.5) {
//            potential->at(s).at(i,j,k) = 0.;
//          } else {
//            potential->at(s).at(i,j,k) = 1.;
//          }	
//	}
//      } 
//    } // end of (i, j, k)-loops
//  } // end of species loop
//}
//// ____________________________________________________________________________
//template <>
//void exp_external_potential_hs(
//		std::vector<BoundaryStructure<3>>& boundary_structures,
//		std::vector<Properties>& species_properties,
//		std::vector<size_t>& affected_species,
//		std::vector<DataFrame<3, double>>* potential) {
//// TODO
//}
//// ____________________________________________________________________________
//template <size_t dim>
//void electrostatic_potential<dim>(
//		PoissonSolver& poisson_solver,
//		DataFrame<dim, double>* potential) {
//}
//template <>
//void electrostatic_potential(
//		PoissonSolver& poisson_solver,
//		DataFrame<3, double>* potential) {
//  // Number of grid points (voxels)
//  std::vector<size_t> grid_counts = potential->size_dim();
//  size_t voxel_count{1};
//  for (auto& gc : grid_counts) { voxel_count *= gc; }
//  // Convert initial guess from DataFrame to the format required by the
//  // PoissonSolver
//  std::vector<double> solution(voxel_count, 0.);
//  size_t index{0};
//  for (size_t i = 0; i < grid_counts.at(0); ++i) {
//    for (size_t j = 0; j < grid_counts.at(1); ++j) {
//      for (size_t k = 0; k < grid_counts.at(2); ++k) {
//        index = i * grid_counts.at(2) * grid_counts.at(1) +
//            j * grid_counts.at(2) + k;
//        solution.at(index) = potential->at(i,j,k);
//      }
//    }
//  }
//  // Solve Poisson equation
//  std::vector<double> rhs(voxel_count, 0.);
//  poisson_solver.solve(rhs, solution);
//  // reconvert to DataFrame format
//  for (size_t i = 0; i < grid_counts.at(0); ++i) {
//    for (size_t j = 0; j < grid_counts.at(1); ++j) {
//      for (size_t k = 0; k < grid_counts.at(2); ++k) {
//        index = i * grid_counts.at(2) * grid_counts.at(1) +
//            j * grid_counts.at(2) + k;
//        potential->at(i, j, k) = solution.at(index);
//      }
//    }
//  }
//}
