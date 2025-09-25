import data_structure as struct
import numpy as np
#import sys

print("main:")
params = ["grid_counts", "max_dev", "potential"]
data = struct.get_reduced_dict(struct.main_data, params)
struct.print_struct(data, ["folder", "num", *params])

potentials = np.array(data, dtype='object')[:,-1]
if len(np.unique(potentials)) != len(potentials):
  print("Warning: Found multiple data files with same potential values in main_data.")

print("planar:")
params = ["grid_counts", "max_dev", "potential"]
data = struct.get_reduced_dict(struct.planar_data, params)
struct.print_struct(data, ["folder", "num", *params], filename='structure_planar.txt')

potentials = np.array(data, dtype='object')[:,-1]
if len(np.unique(potentials)) != len(potentials):
 print("Warning: Found multiple data files with same potential values in planar_data.")
 
#print("gc:")
#params = ["grid_counts", "max_dev", "potential"]
#data = struct.get_reduced_dict(struct.gc_data, params)
#struct.print_struct(data, ["folder", "num", *params])
