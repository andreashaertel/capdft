# SPDX-FileCopyrightText: 2026 Fabienne Dressler <fab.dressler@web.de>
# SPDX-License-Identifier: LGPL-3.0-or-later
# Plot charge distribution profiles of the cylinder's surface charge.

# ______________________________________________________________________________
# Imports
# ______________________________________________________________________________
import matplotlib.pyplot as plt
import numpy as np
import sys
import load_data_local as dat
sys.path.append("/home/fd160/Documents/capdft/examples/plots/src")
import load_data as DAT

# ______________________________________________________________________________
# Definitions
# ______________________________________________________________________________
folders = sys.argv[1:]
electrolyte_default = "full"

# ______________________________________________________________________________
# Plot slice along z direction
# ______________________________________________________________________________
for i,direc in enumerate(folders):
    data = np.loadtxt(direc+f"/surface_charge_{electrolyte_default}.dat")
    system_lengths = dat.get_system_lengths(direc)
    data = dat.cylinder_slice_z(data, system_lengths[0]/2, system_lengths[1]/2)
    if len(data) == 0:
        print(f"missing data at x={system_lengths[0]/2}nm for {direc}")
    else:
        z = data[:,2]
        charge_density = data[:,3]
        plt.plot(z - system_lengths[2]/2, charge_density, "x", label=direc)
if len(folders) > 1:
    plt.legend()
plt.xlabel("$z - z_0$ [nm]")
plt.ylabel("surface charge density [e/nm²]")
plt.title("Slice along z-direction")
plt.grid()
plt.show()

# ______________________________________________________________________________
# Plot slice through the cylinder mantle along a horizontal plane
# ______________________________________________________________________________
for offset in [0,5.5]:
    for i,direc in enumerate(folders):
        data = np.loadtxt(direc+f"/surface_charge_{electrolyte_default}.dat")
        system_lengths = dat.get_system_lengths(direc)
        z = system_lengths[2]/2 + offset
        data, phi = dat.cylinder_slice_phi(data, z, system_lengths[0]/2, system_lengths[1]/2)
        if len(data) == 0:
            print(f"missing data at z={z} for {direc}")
        else:
            charge_density = data[:,3]
            plt.plot(phi, charge_density, "x", label=f"$z - z_0 = {offset}$nm ({direc})")
plt.legend()
plt.xlabel(r"$\phi$")
plt.ylabel("surface charge density [e/nm²]")
plt.title("Slice through cylinder mantle along horizontal planes")
plt.grid()
plt.show()

# ______________________________________________________________________________
# Plot slice over the cylinder cap along x-direction
# ______________________________________________________________________________
for i,direc in enumerate(folders):
    data = np.loadtxt(direc+f"/surface_charge_{electrolyte_default}.dat")
    system_lengths = dat.get_system_lengths(direc)
    z = system_lengths[2]/2 + offset
    data = dat.cylinder_slice_x(data, ymid=system_lengths[1]/2, zmid=system_lengths[2]/2)
    if len(data) == 0:
        print(f"missing data at y={system_lengths[1]/2} for {direc}")
    else:
        x = data[:,0] - system_lengths[0]/2
        y = data[:,1] - system_lengths[1]/2
        charge_density = data[:,3]
        plt.plot(x, charge_density, "x", label=f"{direc}")
if len(folders) > 1:
    plt.legend()
plt.xlabel(r"$x - x_0$ [nm]")
plt.ylabel("surface charge density [e/nm²]")
plt.title("Slice over cylinder cap along x-direction")
plt.grid()
plt.show()
