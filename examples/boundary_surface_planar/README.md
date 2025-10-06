# Examples with planar boundary surfaces

There are three options to model planar boundary surfaces with this code:
1. In one dimension: using PlanarPoissonSolver, FunctionalFMTPlanar etc.
2. In three dimensions: using CartesianPoissonSolver, FunctionalFMTCartesian etc.
3. In three-dimensions: using BoundarySurfacePlanar, CartesianPoissonSolverAny, FunctionalFMTCartesian etc.
The first option is the easiest and fastest one. Thus, it is usually the
preferred choice. However, trying the other options and comparing their
results allows to check the different algorithms' mutual consistency and
correct performance.

This folder contains example code for all three options.
