# SPDX-FileCopyrightText: 2026 Fabienne Dressler <fab.dressler@web.de>
# SPDX-License-Identifier: LGPL-3.0-or-later
# Define methods to extract data from the local parameter and data files.

# ______________________________________________________________________________
# Imports
# ______________________________________________________________________________
import numpy as np
import sys
sys.path.append("/home/fd160/Documents/capdft/examples/plots/src")
import load_data as DAT

# ______________________________________________________________________________
# Definitions
# ______________________________________________________________________________
default_binsize = 0.08 # [nm], valid for most runs (cf. README.md)
electrolyte_default = "full" # alternative: "vacuum" (i.e. ignore ions)

# ______________________________________________________________________________
# System properties
# ______________________________________________________________________________
def get_system_lengths(folder):
    return np.array([DAT.get_param(f"{i}_length", folder) for i in range(3)])

def get_cap_type(folder):
    return DAT.get_param("0_caps", folder).replace(" ","").replace("\n","")

# ______________________________________________________________________________
# Results (scalar values - e.g. averages)
# ______________________________________________________________________________
def get_charges(folder, region : ["cap", "mantle", "total"],
                electrolyte : ["full", "vacuum"] = electrolyte_default,
                density=False):
    # specify where to find the requested information
    label = ""
    if region == "cap":
        label = "# on upper cap:"
    elif region == "mantle":
        label = "# on cylinder mantle:"
    elif region == "total":
        label = "# total surface area:"
    if density:
        in_line_position = -1
        cutoff_unit = -6
    else:
        in_line_position = 1
        cutoff_unit = -1
    # extract value
    string = DAT.get_param(label, 
                           folder+f"/surface_charge_{electrolyte}.dat", 
                           parameter_file=False).split(" ")[in_line_position]
    value = float(string[:cutoff_unit])
    return value

def yukawa_fit_name(version):
    """Specify in which file the results from a certain fit can be found"""
    if version == "simple":
        return "yukawa_fit.dat"
    elif version == "detailed":
        return "yukawa_detailed_fit.dat"
    else:
        raise ValueError("yukawa_fit_name: Unknown version {version}")

def get_yukawa_amplitudes(folder, version="simple"):
    # specify where to find the requested information
    file = f"{folder}/" + yukawa_fit_name(version)
    labels = []
    if version == "simple":
        labels = ["# A: ", "# A_offset: "]
    else:
        number_of_spheres = 6
        labels = [f"# A{i}: " for i in range(number_of_spheres)]
    # extract values
    values = []
    for label in labels:
        string = DAT.get_param(label, file, parameter_file=False)
        value, unit = string.split(" ")
        values.append(float(value))
    return values

def get_yukawa_screening(folder, version="simple"):
    # specify where to find the requested information
    file = f"{folder}/" + yukawa_fit_name(version)
    label = "# l: "
    # extract values
    string = DAT.get_param(label, file, parameter_file=False)
    value, unit = string.split(" ")
    return float(value)

# ______________________________________________________________________________
# Data profiles
# ______________________________________________________________________________
def cylinder_slice_z(dataset, xmid, ymid, acceptance=0.25*default_binsize):
    selection = []
    for data in dataset:
        if abs(data[0] - xmid) < acceptance and data[1] >= ymid:
            selection.append(data)
    return np.array(selection)

def cylinder_slice_phi(dataset, z, xmid, ymid, acceptance=0.25*default_binsize):
    selection = []
    phi = []
    while len(selection) == 0:
        selection = []
        for data in dataset:
            if abs(data[2] - z) < acceptance:
                selection.append(data)
                phi.append(np.arctan((data[1]-ymid)/(data[0]-xmid)))
        acceptance *= 2 # repeat with higher acceptance if no data were found
    return np.array(selection), np.array(phi)

def cylinder_slice_x(dataset, ymid, zmid, acceptance=0.25*default_binsize):
    selection = []
    for data in dataset:
        if abs(data[1] - ymid) < acceptance and data[2] > zmid:
        #if data[2] > zmid:
            selection.append(data)
    return np.array(selection)

# ______________________________________________________________________________
# Miscellaneous
# ______________________________________________________________________________
#def average_x_y(dataset, column):
#    """Assuming that dataset is ordered by z-values (axis 2)"""
#    zvalues = []
#    result = [0]
#    z = dataset[0,2]
#    for data in dataset:
#        if data[2] > z:
#            z = data[2]
#            result.append(data[column])
#
#    return np.array(selection)

#class Cylinder():
#    def __init__(self, xmid=None, ymid=None, zmid=None):
#        self.xmid = xmid
#        self.ymid = ymid
#        self.zmid = zmid
#    def load(self, folder, index=0):
#        self.xmid = DAT.get_param(f"{index}_x_mid", folder, default = self.xmid)
#        self.ymid = DAT.get_param(f"{index}_y_mid", folder, default = self.ymid)
#        self.zmid = DAT.get_param(f"{index}_z_mid", folder, default = self.zmid)
#        self.radius = DAT.get_param(f"{index}_radius", folder)
#        return
#    def __repr__(self):
#        return (f"Cylinder located at ({self.xmid},{self.ymid},{self.zmid}) "
#                f"with radius {self.radius}")


