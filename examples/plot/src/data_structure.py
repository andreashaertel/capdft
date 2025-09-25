# File documenting the data structure in this folder.
#______________________________________________________________________________
# Includes
import sys
import numpy as np
#sys.path.insert(1, '~/Documents/capdft/examples/plot/')
#sys.path.insert(1, '~/Uebungen/capdft-master/examples/plot/')
import load_data as dat

#______________________________________________________________________________
# Definitions
MAX_NUM = 100 # maximum number of runs per folder
# folder structure
struct_direc = "../../undulated/out/"
#planar_direc = "~/Uebungen/capdft-master/examples/test_potentials_planar/out/"
planar_direc = "../../planar/out/"
main_data = [struct_direc+f"main/run0{num}" for num in range(1,9)]
example = (struct_direc+"main/run04", 1) # 0.1V, reference for the following
gc_data = [struct_direc+f"tests/corr_{name}" for name in ["gc_test1","gc_test2","higher_gc"]]
# planar data from external directory
planar_data = [planar_direc + "main/run3"]
# planar_gc_data = [planar_direc + "archive/" + folder for folder in ["run1","run2","test_dist"]]
def is_planar(folder): # identify corresponding directory
  return folder[:len(planar_direc)] == planar_direc
# file names
def densities(folder, num):
  if is_planar(folder): return folder + f"/run{num}_1d_profiles.dat"
  else: return folder + f"/run{num}_3d_profiles.dat"
def potential(folder, num):
  if is_planar(folder): return folder + f"/run{num}_1d_total_ES.dat"
  else: return folder + f"/run{num}_total_ES.dat"
def electrode_charge(folder, num, method):
  if method in ["Poisson", "integral", "diff. (order 1)"]:
    return folder + f"/run{num}_torder1num_electrode_charge.dat"
  elif method == "diff. (order 2)":
    return folder + f"/run{num}_torder2num_electrode_charge.dat"
  elif method == "diff. (order 1 & ana.)":
    return folder + f"/run{num}_torder1ana_electrode_charge.dat"
  elif method == "diff. (order 2 & ana.)":
    return folder + f"/run{num}_torder2ana_electrode_charge.dat"
  else: print(f"invalid method '{method}'")

#______________________________________________________________________________
# Functions
# check which iterations have converged
def is_converged(folder, num):
  max_dev = dat.check_iter(folder, num)
  if eval(max_dev) > 1e-6 or max_dev == "nan":
    return False
  elif folder == struct_direc + "main/run08" and num > 5:
  #I'm excluding these from the analysis because the iterator did some strange
  #things there: look at the iterator.dat-files to see that in the last few
  #calculations in ../undulated/out/main/run08, all yield exactly the same
  #results since the iteration aborted immediately.
    return False
  else:
    return True

def are_converged(folder):
  converged = []
  for num in range(0,MAX_NUM):
    try: conv = is_converged(folder, num)
    except: conv = False
    if conv: converged.append(num)
  return converged

# get an overview of the data parameters
def get_dict(folders, params):
  struct = []
  for folder in folders:
    entry = [folder]
    for param in params:
      if param == "conv": value = are_converged(folder)
      else: value = dat.get_param(param, folder)
      entry.append(value)
    struct.append(entry)
  return np.array(struct, dtype='object')

def get_reduced_dict(folders, params):
  struct = []
  for folder in folders:
    converged = are_converged(folder)
    for num in converged:
      entry = [folder, num]
      for param in params:
        if param == "max_dev":
          value = dat.check_iter(folder, num)
        else:
          value = dat.get_param(param, folder)
          if param == "potential": # extract num-specific value
            try: value = value[num]
            except: pass
        entry.append(value)
      struct.append(entry)
  return np.array(struct, dtype='object')

# print dictionary
def print_struct(struct, labels, filename='structure.txt', header=None):
  f = open(filename, 'w')
  f.write("# ")
  for label in labels:
    f.write(f"[{label}] ")
  f.write("\n")
  for entry in struct:
    for data in entry:
      print(data, file=f, end=' ')
    f.write("\n")
  print("see info in " + filename)

def potential_data(potentials, base=main_data):
  struct = get_reduced_dict(base, ["potential"])
  corr_files = []
  if type(potentials) == float:
    potentials = [potentials]
  for corr_pot in potentials:
    print("search ", corr_pot)
    found = False
    for f,n,p in struct:
      print(p)
      #if p == corr_pot:
      if abs(p-corr_pot) < 1e-6:
        found = True
        corr_files.append([f,n])
        break
    if not found:
      print("could not find potential ", corr_pot)
      exit(1)
  return np.array(corr_files, dtype='object')

import data_structure as STR
import electrode_charge_density as ECD
import load_data as DAT

def slice_1d(meshes, direc, pos):
  slices = []
  for mesh in meshes:
    mesh = np.array(mesh)
#    print(mesh.shape)
    if direc == 0:
      ind = int(pos*mesh.shape[1])
      slices.append(mesh[:,ind])
    elif direc == 1:
      ind = int(pos*mesh.shape[0])
      slices.append(mesh[ind,:])
    else: print("slice: invalid direc")
  return slices

class data_set:
  def __init__(self, folders, nums, potentials):
    self.files = [(f,n) for f, n in zip(folders, nums)]
    self.potentials = potentials
    self.len = len(folders)
  def get_gc(self):
    self.gc = []
    for f,n in self.files:
      self.gc.append(DAT.get_param("grid_counts", f))
  def get_densities_xz(self, species=[0,1]):
    self.x_mesh = []
    self.z_mesh = []
    self.densities_xz = [[] for _ in species]
    for f,n in self.files:
      filename = STR.densities(f,n)
      cols = [nu + 3 for nu in species]
      x,z,rho = DAT.load_slice_plane(filename, 0,2,cols)
      self.x_mesh.append(x)
      self.z_mesh.append(z)
      for nu in species:
        self.densities_xz[nu].append(rho[nu])
  def get_densities_z(self, species=[0,1]):
    self.z_vals = []
    self.densities_z_vals = [[] for _ in species]
    for f,n in self.files:
      filename = STR.densities(f,n)
      cols = [nu + 1 for nu in species]
      z,rho = DAT.load_1d(filename, cols)
      self.z_vals.append(z)
      for nu in species:
        self.densities_z_vals[nu].append(rho[nu])
  def get_total_ES_xz(self):
    self.x_mesh = []
    self.z_mesh = []
    self.total_ES_xz = []
    for f,n in self.files:
      filename = STR.potential(f,n)
      x,z,phi = DAT.load_slice_plane(filename, 0,2,[3])
      self.x_mesh.append(x)
      self.z_mesh.append(z)
      self.total_ES_xz.append(phi[0])
  def get_total_ES_z(self):
    self.z_vals = []
    self.total_ES_z_vals = []
    for f,n in self.files:
      filename = STR.potential(f,n)
      z,phi = DAT.load_1d(filename, [1])
      self.z_vals.append(z)
      self.total_ES_z_vals.append(phi[0])
  def get_local_ECD(self, methods):
    self.local_ECD_methods = methods
    self.x_vals = []
    self.ECD_left = [[] for _ in methods]
    self.ECD_right = [[] for _ in methods]
    for f,n in self.files:
      for m in range(len(methods)):
        x,sigma = ECD.get_ECD_profile_x(f, n, methods[m])
        self.ECD_left[m].append(sigma[0])
        self.ECD_right[m].append(sigma[1])
      self.x_vals.append(x)
  def get_mean_ECD(self, methods):
    self.mECD_methods = methods
    self.mECD_left = [[] for _ in methods]
    self.mECD_right = [[] for _ in methods]
    for f,n in self.files:
      for m in range(len(methods)):
        sigma = ECD.get_mean_ECD(f,n,methods[m])
        self.mECD_left[m].append(sigma[0])
        self.mECD_right[m].append(sigma[1])
  def densities_z(self, nu, pos):
    try: return self.densities_z_vals[nu]
    except: return slice_1d(self.densities_xz[nu], 1, pos)
  def total_ES_z(self, pos):
    try: return self.total_ES_z_vals
    except: return slice_1d(self.total_ES_xz, 1, pos)
  def z(self, d):
      try: return self.z_vals[d]
      except: return np.unique(self.z_mesh[d])
  def x(self, d):
      try: return self.x_vals[d]
      except: return np.unique(self.x_mesh[d])
  def local_ECD(self, side):
    if side == "left": return self.ECD_left
    if side == "right": return self.ECD_right
  def mean_ECD(self, side):
    if side == "left": return self.mECD_left
    if side == "right": return self.mECD_right
