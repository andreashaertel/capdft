// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file cartesian_poisson_solver_any.cpp
 *  \brief Source file for the CartesianPoissonSolverAny class.
 *
 *  The file contains the definitions of the CartesianPoissonSolverAny class.
 */
#include "cartesian_poisson_solver_any.hpp"  // NOLINT
#include "boundary_surface.hpp"
#include "boundaries.hpp"
#include "system.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include <stdexcept>
#include "stl_algorithms.hpp"
// _____________________________________________________________________________
CartesianPoissonSolverAny::CartesianPoissonSolverAny() {
}
// _____________________________________________________________________________
CartesianPoissonSolverAny::CartesianPoissonSolverAny(System<3>& system,
    BoundarySurface<3>* boundary_surface) : SparseMatrix(
    system.grid_counts.at(0) * system.grid_counts.at(1) * system.grid_counts.at(2)) {
  this->system = system;
  this->bin_count = system.grid_counts;
  this->bin_size = system.bin_sizes;
  this->periodic_boundaries = system.periodic_boundaries;
  // Calculate the square of the bin sizes
  bin_size_squared.resize(bin_size.size());
  for (size_t i = 0; i < bin_size.size(); ++i) {
    bin_size_squared.at(i) = bin_size.at(i) * bin_size.at(i);
  }
  set_laplacian();
  set_boundary_values(boundary_surface);
}
// _____________________________________________________________________________
CartesianPoissonSolverAny::CartesianPoissonSolverAny(System<3>& system,
    Boundaries<3>* boundaries)
  : SparseMatrix(
    system.grid_counts.at(0) * system.grid_counts.at(1) * system.grid_counts.at(2)) {
  this->system = system;
  this->bin_count = system.grid_counts;
  this->bin_size = system.bin_sizes;
  this->periodic_boundaries = system.periodic_boundaries;
  // Calculate the square of the bin sizes
  bin_size_squared.resize(bin_size.size());
  for (size_t i = 0; i < bin_size.size(); ++i) {
    bin_size_squared.at(i) = bin_size.at(i) * bin_size.at(i);
  }
  set_laplacian();
  set_boundary_values(boundaries);
}
// _____________________________________________________________________________
CartesianPoissonSolverAny::~CartesianPoissonSolverAny() {
}
// _____________________________________________________________________________
void CartesianPoissonSolverAny::solve(
    std::vector<double>& rhs, std::vector<double>& solution) {
  // Termination condition: maximum allowed residual norm
  double accuracy = pow(std::numeric_limits<double>::epsilon(), 1. / 3.);
  // The algorithm becomes inefficient in terms of memory if too many (>50)
  // iterations are performed.
  size_t max_iterations{100};
  size_t iteration_count{0};
  double deviation{0.};
  // Dummy right-hand side with certain elements removed
  add_boundary_values(rhs);
  std::vector<double> rhs_dummy(rhs);
  remove_boundary_values(rhs);
  stl_algorithm::erase(boundary_positions, rhs_dummy);
  // Dummy solution with reduced size
  std::vector<double> solution_dummy(solution);
  stl_algorithm::erase(boundary_positions, solution_dummy);
  // Solver loop
  while (!SparseMatrix::solve(
      rhs_dummy, solution_dummy, max_iterations, accuracy, &deviation)) {
    iteration_count += max_iterations;
    std::cout << "CartesianPoissonSolverAny::solve(): \"";
    std::cout << "iterations: " << iteration_count;
    std::cout << "; deviation: " << deviation << " > " << accuracy << "\"";
    std::cout << std::endl << "\033[A\033[K";
  }
  // Add the known potential values at the right place
  for (size_t i = 0; i < boundary_points.size(); ++i) {
    solution_dummy.insert(
        solution_dummy.begin() + boundary_points.at(i).first,
        boundary_points.at(i).second);
  }
  solution = solution_dummy;
}
// _____________________________________________________________________________
void CartesianPoissonSolverAny::solve(DataFrame<3, double>& rhs,
		DataFrame<3, double>& solution) {
  size_t voxel_count = bin_count.at(0) * bin_count.at(1) * bin_count.at(2);
  std::vector<double> rhs_vector(voxel_count);
  std::vector<double> solution_vector(voxel_count);
  size_t index;
  for (size_t i = 0; i < bin_count.at(0); i++) {
    for (size_t j = 0; j < bin_count.at(1); j++) {
      for (size_t k = 0; k < bin_count.at(2); k++) {
	coordinates_to_index(i, j, k, &index);
	rhs_vector.at(index) = rhs.at(i, j, k);
      }
    }
  }
  solve(rhs_vector, solution_vector);
  for (size_t i = 0; i < bin_count.at(0); i++) {
    for (size_t j = 0; j < bin_count.at(1); j++) {
      for (size_t k = 0; k < bin_count.at(2); k++) {
	coordinates_to_index(i, j, k, &index);
	solution.at(i, j, k) = solution_vector.at(index);
      }
    }
  }
}
// _____________________________________________________________________________
void CartesianPoissonSolverAny::set_laplacian() {
  set_bare_laplacian();
  set_boundary_conditions();
}
// _____________________________________________________________________________
void CartesianPoissonSolverAny::set_bare_laplacian() {
  size_t row{0}, col{0};
  double diagonal_value{0.};
  std::vector<double> offdiag_values(0);
  // Calculate all possible matrix values
  diagonal_value = -2. / bin_size_squared.at(0) - 2. / bin_size_squared.at(1) -
      2. / bin_size_squared.at(2);
  for (auto& bss : bin_size_squared) {
    offdiag_values.push_back(1. / bss);
  }
  // Convert the 3D cartesian laplacian into a matrix
  for (size_t i = 0; i < bin_count.at(0); ++i) {
    for (size_t j = 0; j < bin_count.at(1); ++j) {
      for (size_t k = 0; k < bin_count.at(2); ++k) {
        coordinates_to_index(i, j, k, &row);
        SparseMatrix::set(row, row, diagonal_value);
        // The six (spacially) adjacent voxels might contribute to the
        // laplacian. If the voxel does not exist, it will be ignored and
        // handled as boundary condition later.
        if(coordinates_to_index(i - 1, j, k, &col)) {
          SparseMatrix::set(row, col, offdiag_values.at(0));
        }
        if(coordinates_to_index(i + 1, j, k, &col)) {
          SparseMatrix::set(row, col, offdiag_values.at(0));
        }
        if(coordinates_to_index(i, j - 1, k, &col)) {
          SparseMatrix::set(row, col, offdiag_values.at(1));
        }
        if(coordinates_to_index(i, j + 1, k, &col)) {
          SparseMatrix::set(row, col, offdiag_values.at(1));
        }
        if(coordinates_to_index(i, j, k - 1, &col)) {
          SparseMatrix::set(row, col, offdiag_values.at(2));
        }
        if(coordinates_to_index(i, j, k + 1, &col)) {
          SparseMatrix::set(row, col, offdiag_values.at(2));
        }
      }
    }
  }
}
// _____________________________________________________________________________
void CartesianPoissonSolverAny::set_boundary_conditions() {
  size_t row{0}, col{0};
  if (periodic_boundaries.at(0)) {  // if periodic
    for (size_t i = 0; i < bin_count.at(1); ++i) {
      for (size_t j = 0; j < bin_count.at(2); ++j) {
        coordinates_to_index(0, i, j, &row);  // one face of cube
        coordinates_to_index(bin_count.at(0) - 1, i, j, &col);  // opposing face
        // PBCs are symmetric:
        // bin on one side is linked to the other side and vice versa
        SparseMatrix::set(row, col, 1. / bin_size_squared.at(0));
        SparseMatrix::set(col, row, 1. / bin_size_squared.at(0));
      }
    }
  } else {  // if not periodic
  }
  if (periodic_boundaries.at(1)) {
    for (size_t i = 0; i < bin_count.at(0); ++i) {
      for (size_t j = 0; j < bin_count.at(2); ++j) {
        coordinates_to_index(i, 0, j, &row);
        coordinates_to_index(i, bin_count.at(1) - 1, j, &col);
        // PBCs are symmetric:
        // bin on one side is linked to the other side and vice versa
        SparseMatrix::set(row, col, 1. / bin_size_squared.at(1));
        SparseMatrix::set(col, row, 1. / bin_size_squared.at(1));
      }
    }
  } else {
  }
  if (periodic_boundaries.at(2)) {
    for (size_t i = 0; i < bin_count.at(0); ++i) {
      for (size_t j = 0; j < bin_count.at(1); ++j) {
        coordinates_to_index(i, j, 0, &row);
        coordinates_to_index(i, j, bin_count.at(2) - 1, &col);
        // PBCs are symmetric:
        // bin on one side is linked to the other side and vice versa
        SparseMatrix::set(row, col, 1. / bin_size_squared.at(2));
        SparseMatrix::set(col, row, 1. / bin_size_squared.at(2));
      }
    }
  } else {
  }
}
// _____________________________________________________________________________
bool CartesianPoissonSolverAny::coordinates_to_index(
    size_t i, size_t j, size_t k, size_t* index) {
  // i, j, k cannot be smaller 0, instead they jump to the largest value.
  if (i >= bin_count.at(0) || j >= bin_count.at(1) || k >= bin_count.at(2)) {
    return false;
  }
  *index = k + bin_count.at(2) * j + bin_count.at(2) * bin_count.at(1) * i;
  return true;
}
// _____________________________________________________________________________
bool CartesianPoissonSolverAny::coordinates_to_index(
    std::vector<size_t> pos, size_t* index) {
  return coordinates_to_index(pos.at(0), pos.at(1), pos.at(2), index);
}
// _____________________________________________________________________________
void CartesianPoissonSolverAny::coordinates_to_position(
    std::vector<size_t> pos, std::vector<double>* position) {
  *position = system.coordinates_to_position(pos);
}
// _____________________________________________________________________________
void CartesianPoissonSolverAny::set_boundary_values(
		Boundaries<3>* boundaries) {
  // Reset global variables
  boundary_points.clear();
  boundary_positions.clear();
  rhs_addition.clear();
  rhs_addition.resize(bin_count.at(0) * bin_count.at(1) * bin_count.at(2), 0.);
  // Calculate modifications (cols/boundary_positions and rhs_addition) for all
  // boundary objects
  for (auto it = boundaries->begin(); it != boundaries->end(); it++) {
    if ((*it)->is_electrostatic()) {
      calc_boundary_values(*it);
    }
  }
  // Sort rows/boundary_points by index
  std::sort(boundary_points.begin(), boundary_points.end());
  std::sort(boundary_positions.begin(), boundary_positions.end());
  // Remove dublicates
  boundary_positions.erase(std::unique(boundary_positions.begin(), boundary_positions.end()),
		  boundary_positions.end());
  boundary_points.erase(std::unique(boundary_points.begin(), boundary_points.end()),
		  boundary_points.end());
  std::cout << "CartesianPoissonSolverAny: " << boundary_points.size()
	  << " surface points\n";
  // What about colliding obj. with different boundary values?
  // There will be more dublicates among the boundary positions than among the
  // position-value pairs.
  if (boundary_positions.size() != boundary_points.size()) {
    std::cerr << "CartesianPoissonSolverAny::set_boundary_values: Error: "
	    "boundary values not well defined - Surfaces with different "
	    "boundary values overlap at " <<
	    boundary_points.size() - boundary_positions.size() << " positions.\n";
  }
  // Remove rows and columns corresponding to any points outside the boundary
  SparseMatrix::remove_columns(boundary_positions);
  SparseMatrix::transpose();
  SparseMatrix::remove_columns(boundary_positions);
  SparseMatrix::transpose();
}
// _____________________________________________________________________________
void CartesianPoissonSolverAny::set_boundary_values(
		BoundarySurface<3>* boundary_surface) {
  // Reset global variables
  boundary_points.clear();
  boundary_positions.clear();
  rhs_addition.clear();
  rhs_addition.resize(bin_count.at(0) * bin_count.at(1) * bin_count.at(2), 0.);
  // Calculate modifications (cols/boundary_positions and rhs_addition) for all
  // boundary objects
  if (boundary_surface->is_electrostatic()) {
    calc_boundary_values(boundary_surface);
  }
  // Remove rows and columns corresponding to any points outside the boundary
  SparseMatrix::remove_columns(boundary_positions);
  SparseMatrix::transpose();
  SparseMatrix::remove_columns(boundary_positions);
  SparseMatrix::transpose();
}
//// _____________________________________________________________________________
//void CartesianPoissonSolverAny::calc_boundary_values(
//		BoundarySurface<3>* boundary_surface) {
//  // Iterate over all grid points outside the boundaries
//  std::vector<double> position = {0., 0., 0.};
//  std::vector<size_t> coordinates = {0, 0, 0};
//  size_t index = 0;
//  double value = boundary_surface->get_boundary_value();
//  for (size_t i = 0; i < bin_count.at(0); ++i) {
//    position.at(0) = static_cast<double>(i) * bin_size.at(0);
//    for (size_t j = 0; j < bin_count.at(1); ++j) {
//      position.at(1) = static_cast<double>(j) * bin_size.at(1);
//      for (size_t k = 0; k < bin_count.at(2); ++k) {
//        position.at(2) = static_cast<double>(k) * bin_size.at(2);
//	bool is_in_boundaries = boundary_surface->is_within_boundary(position);
//	if (!is_in_boundaries) {
//	  coordinates = {i, j, k};
//	  coordinates_to_index(coordinates, &index);
//	  // Modify rhs and (off-)diagonal values of adjacent grid points, 
//	  // if they lie within boundaries
//	  modify_rhs(boundary_surface, coordinates);
//	  // Remove corresponding row and column. The removed indices and
//	  // values are stored so that they can be filled back in to the final
//	  // solution vector.
//	  boundary_points.push_back(std::make_pair(index, value));
//	  boundary_positions.push_back(index);
//	}
//      }
//    }
//  }
//}
// _____________________________________________________________________________
void CartesianPoissonSolverAny::calc_boundary_values(
		BoundarySurface<3>* boundary_surface) {
  // Iterate over all grid points outside the boundaries
  std::vector<double> position = {0., 0., 0.};
  std::vector<size_t> coordinates = {0, 0, 0};
  size_t index = 0;
  double value = boundary_surface->get_boundary_value();
  for (size_t i = 0; i < bin_count.at(0); ++i) {
    position.at(0) = static_cast<double>(i) * bin_size.at(0);
    for (size_t j = 0; j < bin_count.at(1); ++j) {
      position.at(1) = static_cast<double>(j) * bin_size.at(1);
      for (size_t k = 0; k < bin_count.at(2); ++k) {
        position.at(2) = static_cast<double>(k) * bin_size.at(2);
	coordinates = {i, j, k};
	bool is_in_boundaries = boundary_surface->is_within_boundary(position);
	if (!is_in_boundaries) {
	  // Remove corresponding row and column. The removed column indices and
	  // values are stored so that they can be filled back in to the final
	  // solution vector.
	  coordinates_to_index(coordinates, &index);
	  boundary_points.push_back(std::make_pair(index, value));
	  boundary_positions.push_back(index);
      } else {
	  // Modify rhs and matrix elements associated with this grid point
	  // depending on its neighbors
	  modify_matrix_elements(boundary_surface, coordinates);
	}
      }
    }
  }
}
// _____________________________________________________________________________
void CartesianPoissonSolverAny::modify_matrix_elements(
		BoundarySurface<3>* boundary_surface,
		std::vector<size_t>& coordinates) {
  double value = boundary_surface->get_boundary_value();
  std::vector<double> position = system.coordinates_to_position(coordinates);
  size_t index;
  size_t index_neighbor;
  std::vector<size_t> coordinates_neighbor(3);
  // Check whether the given coordinate is in direct vicinity (closer than bin
  // size) to the surface
  double distance_left;
  double distance_right;
  double rhs;
  for (size_t dir = 0; dir < 3; dir++) {
    distance_left = boundary_surface->distance_directed(position, dir, false);
    distance_right = boundary_surface->distance_directed(position, dir, true);
    /* We are only interested in whether the distance is smaller than the bin
     * size. If it is larger than that, we set it to -1, meaning there is no
     * surface in vicinity.
     */
    if (distance_left > bin_size.at(dir)) {
      distance_left = -1.;
    }
    if (distance_right > bin_size.at(dir)) {
      distance_right = -1.;
    }
    // Initialize rhs addition
    rhs = 0.;
    if (distance_left < 0. && distance_right < 0.) {
      // No surface nearby: do nothing, check other directions
    } else {
      // Surface adjacent to this point: modify matrix elements in this row
      coordinates_to_index(coordinates, &index);
      coordinates_neighbor = coordinates;
      if (distance_left > 0.) {
	// Modify rhs
	rhs += 1. / distance_left;
      } else {
	// Set default value for distance to neighbor (for calculation of rhs &
	// diagonal matrix element (index, index)); do not change anything else
	distance_left = bin_size.at(dir);
      }
      if (distance_right > 0.) {
	// Modify rhs
	rhs += 1. / distance_right;
      } else {
	// Set default value for distance to neighbor
	distance_right = bin_size.at(dir);
      }
      /** Modify rhs and matrix elements in the row corresponding to this point
       */
      /* Right hand side modification:
       * TODO: generalize for overlapping surfaces - this here would count it twice
       */
      rhs_addition.at(index) -=
		  value * 2. / (distance_left + distance_right) * rhs;
      /* Diagonal matrix element (index, index):
       */
      double old_value = -2. / bin_size_squared.at(dir);
      double new_value = -2. / (distance_left * distance_right);
      SparseMatrix::set(index, index,
			SparseMatrix::get(index, index) - old_value + new_value);
      /* Off-diagonal matrix elements (index, index_neighbor):
       * Zero if there's a boundary inbetween the two points, else nonzero.
       */
      try { // left neighbor
        coordinates_neighbor = coordinates;
        coordinates_neighbor.at(dir) = system.increase(coordinates.at(dir), dir, -1);
        coordinates_to_index(coordinates_neighbor, &index_neighbor);
        if (distance_left < bin_size.at(dir)) {
          SparseMatrix::set(index, index_neighbor, 0.);
        } else {
	  SparseMatrix::set(index, index_neighbor,
		2. / (distance_left * (distance_left + distance_right)));
	}
      } catch (const std::out_of_range&) {
        // do nothing if there is no neighbor in this direction
      }
      try { // right neighbor
        coordinates_neighbor = coordinates;
        coordinates_neighbor.at(dir) = system.increase(coordinates.at(dir), dir, +1);
        coordinates_to_index(coordinates_neighbor, &index_neighbor);
        if (distance_right < bin_size.at(dir)) {
          SparseMatrix::set(index, index_neighbor, 0.);
        } else {
	  SparseMatrix::set(index, index_neighbor,
		2. / (distance_right * (distance_left + distance_right)));
	}
      } catch (const std::out_of_range&) {
        // do nothing if there is no neighbor in this direction
      }
    }
  }
}
//// _____________________________________________________________________________
//void CartesianPoissonSolverAny::modify_rhs(
//		BoundarySurface<3>* boundary_surface,
//		std::vector<size_t>& coordinates1) {
//  // index & boundary value at given point (must be a point outside boundaries!)
//  size_t index1;
//  coordinates_to_index(coordinates1, &index1);
//  double value = boundary_surface->get_boundary_value();
//  // indices & coordinates of next & next-to-next neighbors in the grid
//  size_t index2, index3;
//  std::vector<size_t> coordinates2(3);
//  std::vector<size_t> coordinates3(3);
//  std::vector<double> position2(3);
//  // Iterate over the six next neighbors
//  for (size_t dir = 0; dir < 3; dir++) {
//    for (int step : {-1, 1}) {
//      try {
//        // Calculate neighbor coordinates position2
//        coordinates2 = coordinates1;
//        coordinates2.at(dir) = system.increase(coordinates1.at(dir), dir, step);
//        coordinates_to_position(coordinates2, &position2);
//        // Calculate distance to boundary in the direction position2 -> position1
//	// (forward if step was backwards, i.e. if step==-1)
//	bool forward = (step == -1) ? true : false;
//        double distance = boundary_surface->distance_directed(position2, dir,
//				forward);
//	// If this neighbor is also within the boundaries, distance_directed returns
//	// zero. In that case, we don't have to do anything. Else:
//        if (distance > 0.) {
//	  // Modify rhs and matrix elements in the row corresponding to index2
//    	  coordinates_to_index(coordinates2, &index2);
//	  // Assuming there is no directly adjacent boundary in the opposite 
//	  // direction (else different old_value & new_value):
//	  // TODO: generalize for cases where two surfaces are closer than 2*bin_size
//          // Matrix element at (index2, index2):
//	  double old_value = -2. / bin_size_squared.at(dir);
//	  double new_value = -2. / (bin_size.at(dir) * distance);
//	  SparseMatrix::set(index2, index2,
//			SparseMatrix::get(index2, index2) - old_value + new_value);
//	  // Matrix element at (index2, index1) is moved to right hand side:
//	  // - old_value(rhs) ~ density at index2
//	  // - old_value(matrix element) = 1. / bin_sizes_squared.at(dir);
//	  // TODO: generalize for overlapping surfaces - this here would count it twice
//          rhs_addition.at(index2) -=
//		  value * 2. / (distance * (bin_size.at(dir) + distance));
//          // Matrix element at (index2, index3):
//	  try {
//	    // Calculate the coordinates of the next-to-next neighbor index3
//            coordinates3 = coordinates2;
//            coordinates3.at(dir) = system.increase(coordinates2.at(dir), dir, step);
//	    coordinates_to_index(coordinates3, &index3);
//	    // - old_value = 1. / bin_size_squared.at(dir)
//	    SparseMatrix::set(index2, index3,
//		2. / (bin_size.at(dir) * (bin_size.at(dir) + distance)));
//	  } catch (const std::out_of_range&) {
//            // do nothing if there is no next-to-next neighbor in this direction
//	  } catch (...) {
//	    std::cout << "CartesianPoissonSolverAny::modify_rhs: Error\n";
//	    exit(1);
//          }
//        }
//      } catch (const std::out_of_range&) {
//        // do nothing if there is no next neighbor in this direction
//      } catch (...) {
//	std::cout << "CartesianPoissonSolverAny::modify_rhs: Error\n";
//	exit(1);
//      }
//    }
//  }
//}
// _____________________________________________________________________________
void CartesianPoissonSolverAny::add_boundary_values(
    std::vector<double>& rhs) {
  for (size_t i = 0; i < rhs.size(); ++i) { rhs.at(i) += rhs_addition.at(i); }
}
// _____________________________________________________________________________
void CartesianPoissonSolverAny::remove_boundary_values(
    std::vector<double>& rhs) {
  for (size_t i = 0; i < rhs.size(); ++i) { rhs.at(i) -= rhs_addition.at(i); }
}
// _____________________________________________________________________________
// _____________________________________________________________________________
