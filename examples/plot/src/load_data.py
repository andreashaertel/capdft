import numpy as np
import matplotlib.pyplot as plt

def is_in_plane(index, vol_shape, P, plane_pos):
    """Filter out in-plane-data based on line number, so that not all lines have to be read."""
    if vol_shape[0] == -1:
        print("no header line specifying data shape has been found. place it above all data lines.")
        exit()
    if plane_pos is None: plane_pos = int(vol_shape[P]/2) # default: middle plane
    else: plane_pos = int(vol_shape[P] * plane_pos) # rescale from [0,1] to [0, grid count]
    index -= 1
    if P == 0: return int(index / (vol_shape[1]*vol_shape[2])) == plane_pos
    if P == 1: return int((index % (vol_shape[1]*vol_shape[2])) / vol_shape[2]) == plane_pos
    if P == 2: return (index % (vol_shape[1]*vol_shape[2])) % vol_shape[2] == plane_pos

def get_vol_shape(line):
  text = line.split('=')
  if text[0] == '# grid_counts':
   return eval(text[1])
 
def get_param(name, direc):
  with open(direc+"/params.txt", 'r') as file:
    for line in file:
      if line[:len(name)] == name:
        vals = line[len(name):]
        vals = vals.split("#")[0]
        return eval(vals)
  print(f"parameter {name} not found in {direc}/params.txt.")
  return None

def check_iter(direc, num):
  try:
    file = open(direc+f"/run{num}_iterator.dat", 'r')
    for line in file:
      if line[0] == "#" or line=="": pass
      else: maxdev = line.split(' ')[0]
  except: pass #print("check_iter: Error")
  return maxdev

# load data:
def load_3d(file, V):
  X,Y,Z = 0,1,2
  x = []
  y = []
  z = []
  v = [[] for _ in V]
  i = 0
  j = 0
  print("reading data from ", file)
  with open(file, 'r') as f:
      for line in f:
          if line[0] == '#':
              j += 1
              print(line[:-1])
              temp = get_vol_shape(line)
              if temp is not None: vol_shape = temp
          elif line == '\n':
              pass
          else:
              i += 1
              if is_in_plane(i, vol_shape, P, plane_pos):
                  vals = line.split(' ')
                  x.append(float(vals[X]))
                  y.append(float(vals[Y]))
                  z.append(float(vals[Z]))
                  for s in range(len(V)):
                      if len(vals) <= max(V):
                          print(f"missing vals in data line {i}: only found {len(vals)} values.")
                          exit(1)
                      val = float(vals[V[s]])
                      v[s].append(val)
  print(f"found {i} data lines and {j} comments")

  # check the shape:
  xyz_shape = (len(np.unique(x)), len(np.unique(y)), len(np.unique(z)))
  v = np.array(v)
  if vol_shape != xyz_shape:
      print("it seems like something has gone wrong: "
              f"expected vol shape {vol_shape} and found {xyz_shape} do not match.")
      exit(1)
  elif vol_shape[0]*vol_shape[1]*vol_shape[2] != v.shape[1]:
      print("it seems like something has gone wrong: "
              f"vol shape {vol_shape} and found values shape {v.shape} do not match.")
      exit(1)
  # reshape
  x = np.reshape(x, vol_shape)
  y = np.reshape(y, vol_shape)
  z = np.reshape(z, vol_shape)
  z = np.reshape(v, (len(V), *vol_shape))
  return x,y,z,v

def load_slice_plane(file, X, Y, Z, plane_pos=None):
  vol_shape = (-1,-1,-1)
  P = 3 - X - Y
  x = []
  y = []
  z = [[] for _ in Z]
  i = 0
  j = 0
  print("reading data from ", file)
  with open(file, 'r') as f:
      for line in f:
          if line[0] == '#':
              j += 1
              print(line[:-1])
              temp = get_vol_shape(line)
              if temp is not None: vol_shape = temp
          elif line == '\n':
              pass
          else:
              i += 1
              if is_in_plane(i, vol_shape, P, plane_pos):
                  vals = line.split(' ')
                  x.append(float(vals[X]))
                  y.append(float(vals[Y]))
                  for s in range(len(Z)):
                      if len(vals) <= max(Z):
                          print(f"missing vals in data line {i}: only found {len(vals)} values.")
                          exit(1)
                      val = float(vals[Z[s]])
                      z[s].append(val)
  print(f"found {i} data lines and {j} comments")

  # check the shape:
  plane_shape = (vol_shape[X], vol_shape[Y])
  xy_shape = (len(np.unique(x)), len(np.unique(y)))
  z = np.array(z)
  if plane_shape != xy_shape:
      print("it seems like something has gone wrong: "
              f"expected plane shape {plane_shape} and found {xy_shape} do not match.")
      print("len(x) = ", len(x))
      exit(1)
  elif plane_shape[0]*plane_shape[1] != z.shape[1]:
      print("it seems like something has gone wrong: "
              f"plane shape {plane_shape} and found z shape {z.shape} do not match.")
      exit(1)
  # reshape
  x = np.reshape(x, plane_shape)
  y = np.reshape(y, plane_shape)
  z = np.reshape(z, (len(Z), *plane_shape))

  return x,y,z
#______________________________________________________________________________
def load_2d(file, X, Y, Z):
  x = []
  y = []
  z = [[] for _ in Z]
  i = 0
  j = 0
  print("reading data from ", file)
  with open(file, 'r') as f:
      for line in f:
          if line[0] == '#':
              j += 1
              print(line[:-1])
          elif line == '\n':
              pass
          else:
              i += 1
              vals = line.split(' ')
              x.append(float(vals[X]))
              y.append(float(vals[Y]))
              for s in range(len(Z)):
                  if len(vals) <= max(Z):
                      print(f"missing vals in data line {i}: only found {len(vals)} values.")
                      exit(1)
                  val = float(vals[Z[s]])
                  z[s].append(val)
  print(f"found {i} data lines and {j} comments")

  # check the shape:
  xy_shape = (len(np.unique(x)), len(np.unique(y)))
  z = np.array(z)
  if xy_shape[0]*xy_shape[1] != z.shape[1]:
      print("it seems like something has gone wrong: "
              f"found xy shape {xy_shape} and z shape {z.shape} do not match.")
      exit(1)
  # reshape
  x = np.reshape(x, xy_shape)
  y = np.reshape(y, xy_shape)
  z = np.reshape(z, (len(Z), *xy_shape))

  return x,y,z

#______________________________________________________________________________
def load_1d(file, Z):
  x = []
  z = [[] for _ in Z]
  i = 0
  j = 0
  print("reading data from ", file)
  with open(file, 'r') as f:
      for line in f:
          if line[0] == '#':
              j += 1
              print(line[:-1])
          elif line == '\n':
              pass
          else:
              i += 1
              vals = line.split(' ')
              x.append(float(vals[0]))
              for s in range(len(Z)):
                if len(vals)-1 <= max(Z):
                  print("missing vals in data line ", i)
                  exit(1)
                val = float(vals[Z[s]])
                z[s].append(val)
  print(f"found {i} data lines and {j} comments")
  return x,z

markers = "x+.^s>*............"

def eff_range(z, epsilon=1e-4):
  z = np.array(z)
  epsilon *= (z.max() - z.min())
  if abs(z.min()) < epsilon:
    eff_min = z[z > 0].min()
    eff_max = z[z < z.max() - epsilon].max()
  else:
    eff_min = z[z > z.min() + epsilon].min()
    eff_max = z[z < z.max() - epsilon].max()
  print("plot z range: ", eff_min, eff_max) 
  return eff_min, eff_max

def var_range(z, nsigma=2):
  z = np.array(z)
  mean = np.mean(z)
  delta_z = nsigma * np.std(z)
  return max([mean - delta_z, z.min()]), min([mean + delta_z, z.max()])

def heatmap(a,b,c,crange,labels):
  if type(crange) == str and crange == "auto": crange = eff_range([c])
  if crange is None: plt.pcolormesh(a,b,c)
  else: plt.pcolormesh(a,b,c,clim=crange)
  plt.colorbar(label=labels[2])
  plt.xlabel(labels[0])
  plt.xlabel(labels[1])
  
def heatmap_xz(x,z,c,crange=None,clabel=""):
  heatmap(z,x,c,crange,["$z$ [nm]", "$x$ [nm]", clabel])

