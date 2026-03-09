<!--
SPDX-FileCopyrightText: 2026 Fabienne Dressler <fab.dressler@web.de>

SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Examples

This folder contains examples that illustrate possible applications for the
CapDFT library.

## Overview
A few of the examples illustrate examples with different boundary shapes in the
3d Cartesian geometry:  
boundary_surface_cylinder: One charged cylinder embedded in an ionic liquid  
boundary_surface_sine: Ionic liquid between two undulating electrodes  
boundary_surface_sphere: Two charged spheres at fixed distance, embedded in
	an ionic liquid  

A few of them illustrate examples in the planar or spherical geometry:  
planar_functionals_all_tools: Ionic liquid between two planar electrodes  
spherical_functionals_all_tools: Ionic liquid around a fixed test particle  
spherical_functionals_some_tools: Hard-sphere liquid around a test particle  

A few illustrate differences between different implementation choices:  
boundary_surface_planar: Comparison between 1d planar and 3d Cartesian geometry  
	This example is equivalent to planar_functionals_all_tools, but contains
	three different implementations for the same task:  
	- cartesian_3d_any.cpp: BoundarySurfacePlanar, CartesianPoissonSolverAny  
	- cartesian_3d_planar.cpp: CartesianPoissonSolver ------- DEPRECATED?  
	- planar_1d.cpp: PlanarPoissonSolver  
cartesian_functionals_structured_boundaries and
cartesian_functionals_planar_boundaries: Comparison between 3d Poisson solvers  
	These two try to model the same thing, but one uses the
	CartesianPoissonSolverAny, the other just CartesianPoissonSolver.
	------- DEPRECATED?  
