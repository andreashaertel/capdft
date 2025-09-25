# This file contains stuff concerning the electrode charge density (EDC) data.
import numpy as np
import matplotlib.pyplot as plt
import data_structure as struct
import load_data as dat
import potentials as POT
import plot_basics as PLOT

ECD_methods = ["Poisson", "diff. (order 1)", "diff. (order 2)", "diff. (order 1 & ana.)", "diff. (order 2 & ana.)", "integral"]
ECD_flags = ["# mean charge density (Poisson): "] + 4 * ["# mean charge density (diff): "] + ["# mean charge density (integral): "]
local_ECDs = [1,2,3,4] # all diff methods

main_ECD_methods = [1] # diff
all_ECD_methods = [0,1,2,5] # without ana
planar_ECD_methods = [1] # diff
all_planar_ECD_methods = [1] # diff

def get_mean_ECD(folder, num, method):
  #if type(method) == int: method = ECD_methods[method]
  filename = struct.electrode_charge(folder, num, ECD_methods[method])
  flag = ECD_flags[method]
  right, left = (-10,-10)
  with open(filename, 'r') as f:
    for line in f:
      if line[:len(flag)] == flag:
        values = line[len(flag):].split(' ')
        left, right = eval(values[0]), eval(values[1])
  return left, right

def sort(crit, values):
  order = np.argsort(crit) # find correct order for crit
#  print(order)
  for v in range(len(values)): # reorder arrays
    if len(values[v]) != len(order):
      print(f"sort: array lengths {len(values[v])} and {len(order)} don't match")
    values[v] = [values[v][i] for i in order]
  return np.array(values)

def deriv(x,y):
  if len(x) != len(y): print(f"sizes {len(x)} and {len(y)}")
  x,y = sort(x,[x,y])
  deltax=x[1:]-x[0:-1]
  deltay=y[1:]-y[0:-1]
  return (x[1:]+x[0:-1])/2, deltay/deltax

def plot_c_diff(potentials,mean_ECD,labels=None,colors=None, fontsize=PLOT.std_fontsize, annotate=False, annotate_pos=None):
  if labels is None: labels = range(len(mean_ECD))
  for method in range(len(mean_ECD)):
    p,c = deriv(potentials,mean_ECD[method])
    if colors is None: col = f"C{method}"+dat.markers[method]
    else: col = colors[method]
    plt.plot(p, abs(c), col, label=labels[method])
    if annotate:
      #plt.text(p[-1], abs(c)[-1], labels[method])
      if annotate_pos is None: annotate_pos = -1
      plt.text(p[annotate_pos], abs(c)[annotate_pos], labels[method])
  plt.grid()
  plt.xlabel("$\psi^0$ [V]", size=fontsize)
  plt.ylabel(r"$c^{diff}$ [e/nm²/V]", size=fontsize)

def get_ECD_profile_x(folder, num, method):
  if method not in local_ECDs: return None
  filename = struct.electrode_charge(folder, num, ECD_methods[method])
  x,y,sigma = dat.load_2d(filename, 0, 1, [2,3])
  plane_shape = (len(np.unique(x)), len(np.unique(y)))
  y_pos = int(0.5 * plane_shape[1])
  return x[:,y_pos], sigma[:,:,y_pos]
  
def plot_ECD_profile_x(x, sigma, sides, labels=None):
  if type(sides) == int: sides = [sides]
  for side in sides:
    for i in range(len(sigma)):
      if labels is None: label = i
      else: label = labels[i]
      plt.plot(x[i], sigma[i][side], f'C{i}'+dat.markers[i]+'-', label=label)
  plt.xlabel("$x$ [nm]")
  plt.ylabel("$\sigma$ [e/nm²/V]")
 
# def plot_ECD_profiles_pot(potential, xs, sigma):
#   potential, xs, sigma = sort(potential,[potential, xs, sigma])
#   pmesh = np.array([potential for x in xs[0]]).transpose()
#   xs = np.array(xs, dtype=np.ndarray)
#   print([len(x) for x in xs])
#   xs = np.array(xs)
#   pmesh = np.array(pmesh)
#   sigma = np.array(sigma)
#   print([len(p) for p in pmesh])
#   print([len(s) for s in sigma])
#   plt.pcolormesh(pmesh, xs, sigma, crange=(0,2))
#   #dat.heatmap(pmesh, xs, sigma)
#   #, crange=None,
#    #       labels=["$\psi^0$ [V]", "$x$ [nm]", 
#     #          "$\sigma$ [e/nm²/V]"])

# plot

def get_mean_ECDs(folders, nums, methods=main_ECD_methods):
  data_left = [[] for _ in methods]
  data_right = [[] for _ in methods]
  for folder, num in zip(folders, nums):
    for i, method in enumerate(methods):
      left, right = get_mean_ECD(folder, num, method)
      data_left[i].append(left)
      data_right[i].append(right)
  return data_left, data_right

def plot_mean_ECD(potentials, data_left, data_right, methods, labels=None, color=None):
  if len(data_left) != len(methods): print("Error: ", len(data_left), len(methods))
  for i in range(len(methods)):
    if color==None: col = f"C{methods[i]}"
    else: col = color
    if labels==None: label = ECD_methods[methods[i]]
    else: label = labels[i]
    plt.plot(potentials, data_left[i], col+dat.markers[methods[i]], label=label)
    plt.plot(potentials, data_right[i], col+dat.markers[methods[i]])
  plt.grid()
  plt.xlabel(r"$\psi^0$ [V]")
  plt.ylabel(r"$\hat\sigma$ [e/nm²]")
