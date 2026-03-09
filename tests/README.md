<!--
SPDX-FileCopyrightText: 2026 Fabienne Dressler <fab.dressler@web.de>

SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Tests

This folder contains tests for some of the functions and classes implemented in
the source folder ../src.

## Overview
boundary_surface: Check / visualize output from the BoundarySurface subclasses  
data_frames: Demonstrate usage of DataFrame  
poisson_solver_cartesian: Check that CartesianPoissonSolverAny gives sensible
	results by comparing with analytic expectation, with results from 
	CartesianPoissonSolver or by qualitative visual check  
properties: Demonstrate usage of Properties  
surface_charge: Check that functions from ../src/surface_charge_distribution.hpp
	give sensible results by comparing with analytic expectation  
system: Contains two implementations:  
	- check_interpolation.cpp: Check that System::interpolate gives sensible
		results  
	- interpolate_data.cpp: Can be applied in practice to interpolate data
		from a given data file  
