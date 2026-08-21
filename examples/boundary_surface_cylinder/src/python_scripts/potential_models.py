# SPDX-FileCopyrightText: 2026 Fabienne Dressler <fab.dressler@web.de>
# SPDX-License-Identifier: LGPL-3.0-or-later
# Define effective potential models for the electrostatic potential.

# ______________________________________________________________________________
# Imports
import numpy as np

# ______________________________________________________________________________
# Yukawa potential
def yukawa1(positions, amplitude, screening, R):
    """Model a point charge at position zero, inducing a Yukawa-type potential:
        V(r) = amplitude * exp(- r * screening)/r.
    params:
        positions: array of positions r=(x,y,z)
        amplitude: potential strength
        screening: inverse of the screening length
        R: cutoff radius, set potential to constant below this cutoff
    """
    V = np.zeros(len(positions))
    R_sq = R**2
    cutoff = np.exp(- R * screening) / R
    for i, position in enumerate(positions):
        r_sq = sum(position**2)
        if r_sq < R_sq:
            V[i] += cutoff
        else:
            r = np.sqrt(r_sq)
            V[i] += np.exp(- r * screening) / r
    return amplitude * V

def yukawa2(positions, amplitude, screening, half_height, R=None):
    """Model the potential around the cylindrical colloid by two point charges
    located at the cylinder tips, each of which induces a Yukawa-type potential:
        V(r) = amplitude * (exp(- r1 * screening)/r1 + exp(-r2 * screening)/r2)
        with r1/2 = |r - (0, 0, +/- half_height)|.
    Set it to a constant value for r1/2 < R = cylinder radius.
    params:
        positions: array of positions r=(x,y,z), where r=(0,0,0) is the cylinder
            center
        amplitude: potential strength
        screening: inverse of the screening length
        half_height: parametrizes cylinder tip positions (0, 0, +/-half_height)
        R: cutoff radius - default is half_height / 10
    """
    V = np.zeros(len(positions))
    R = R or (half_height / 10)
    R_sq = R**2
    cutoff = np.exp(- R * screening) / R
    for i, position in enumerate(positions):
        distance_xy_sq = position[0]**2 + position[1]**2
        distance_z_sq = ((position[2] - half_height)**2,
                           (position[2] + half_height)**2)
        for dz_sq in distance_z_sq:
            r_sq = distance_xy_sq + dz_sq
            if r_sq < R_sq:
                V[i] += cutoff
            else:
                r = np.sqrt(r_sq)
                V[i] += np.exp(- r * screening) / r
    return amplitude * V

def yukawa11(positions, amplitude, screening, half_height, R=None):
    """Model the potential around the cylindrical colloid by eleven point charges,
    each of which induces a Yukawa-type potential:
        V(r) = amplitude * sum_i exp(- ri * screening) / ri
        with ri = |r - (0, 0, -half_height + i * 2R)|, i=0,...,10.
    Set it to a constant value for ri < R = cylinder radius.
    params:
        positions: array of positions r=(x,y,z), where r=(0,0,0) is the cylinder
            center
        amplitude: potential strength
        screening: inverse of the screening length
        half_height: parametrizes cylinder tip positions (0, 0, +/-half_height)
        R: cutoff radius - default is half_height / 10
    """
    V = np.zeros(len(positions))
    R = R or (half_height / 10)
    R_sq = R**2
    cutoff = np.exp(- R * screening) / R
    for i, position in enumerate(positions):
        distance_xy_sq = position[0]**2 + position[1]**2
        distance_z = position[2] - half_height
        for j in range(11):
            r_sq = distance_xy_sq + distance_z**2
            if r_sq < R_sq:
                V[i] += cutoff
            else:
                r = np.sqrt(r_sq)
                V[i] += np.exp(- r * screening) / r
            distance_z += 2*R
    return amplitude * V

# ______________________________________________________________________________
# Coulomb potential
def coulomb2(positions, amplitude, half_height):
    """Model the potential around the cylindrical colloid by two point charges
    located at the cylinder tips, each of which induces a Coulomb potential:
        V(r) = amplitude * (1 / r1 + 1/r2)
        with r1/2 = |r - (0,0,+/-half_height)|.
    Set it to a constant value for ri < R = cylinder radius.
    params:
        positions: array of positions r=(x,y,z), where r=(0,0,0) is the cylinder
            center
        amplitude: potential strength
        half_height: parametrizes cylinder tip positions (0, 0, +/-half_height)
        R: cutoff radius - default is half_height / 10
    """
    V = np.zeros(len(positions))
    R = R or (half_height / 10)
    R_sq = R**2
    cutoff = 1 / R
    for i, position in enumerate(positions):
        distance_xy_sq = position[0]**2 + position[1]**2
        distance_z_sq = ((position[2] - half_height)**2,
                         (position[2] + half_height)**2)
        for dz_sq in distance_z_sq:
            r_sq = distance_xy_sq + dz_sq
            if r_sq < R_sq:
                V[i] += cutoff
            else:
                r = np.sqrt(r_sq)
                V[i] += 1 / r
    return amplitude * V

def coulomb11(positions, amplitude, half_height):
    """Model the potential around the cylindrical colloid by eleven point charges,
    each of which induces a Coulomb potential:
        V(r) = amplitude * sum_i 1 / ri
        with ri = |r - (0, 0, -half_height + i * 2R)|, i=0,...,10.
    Set it to a constant value for ri < R = cylinder radius.
    params:
        positions: array of positions r=(x,y,z), where r=(0,0,0) is the cylinder
            center
        amplitude: potential strength
        half_height: parametrizes cylinder tip positions (0, 0, +/-half_height)
    """
    V = np.zeros(len(positions))
    R = half_height / 10
    R_sq = R**2
    cutoff = 1 / R
    for i, position in enumerate(positions):
        distance_xy_sq = position[0]**2 + position[1]**2
        distance_z = position[2] - half_height
        for j in range(11):
            r_sq = distance_xy_sq + distance_z**2
            if r_sq < R_sq:
                V[i] += cutoff
            else:
                r = np.sqrt(r_sq)
                V[i] += 1 / r
            distance_z += 2*R
    return amplitude * V

