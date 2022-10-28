// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file flags.hpp
 *  
 *  This header file contains all flag definitions, such that they can be used
 *  throughout the library.
 */
/** \brief Flag definition for boundary conditions of the PlanarPoissonSolver
 */
enum BoundaryFlag {
  NEUMANN_NEUMANN = 0,
  DIRICHLET_DIRICHLET = 1,
  NEUMANN_DIRICHLET = 2,
  DIRICHLET_NEUMANN = 3
};
