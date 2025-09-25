import numpy as np
import matplotlib.pyplot as plt
import sys
import load_data as dat
add = -100

def load_charge(filename):
#  flags = ["# total charge (diff): ", "# total charge (Poisson): ",
#    "# total charge (integral): "]
  flags = ["# mean charge density (diff): ", "# mean charge density (Poisson): ",
    "# mean charge density (integral): "]
  left = [add] * len(flags)
  right = [add] * len(flags)
  with open(filename, 'r') as f:
    for line in f:
      for i in range(len(flags)):
        flag = flags[i]
        if line[:len(flag)] == flag:
          values = line[len(flag):].split(' ')
          left[i] = float(values[0])
          right[i] = float(values[1])
  return left, right

# default_direcs = ["run21//.","run23//.","run24//.","run25//.","run26//.","run27//.","gc_test1//.","gc_test2//.","higher_gc//.","new_run23//."]
default_direcs = ["corr_" + folder for folder in ["run20//.","run21//.","run23//.","run24//.","run27//.","run28//.","run29//.","run30//.","run31//.","run33//.","run33b//.","run33c//.","run34//.","gc_test1//.","gc_test2//.","higher_gc//.","test_planar//."]]
planar_direcs = ["run1//.","run3//.","test_dist//."]
planar = False
highlighting = False
highlight = 0.1
#colors = ['c', 'k', 'r']

direcs = default_direcs
    #print(f"Usage: {sys.argv[0]} <directory0/>;<directory1>;... <options>\n"
    #       "or {sys.argv[0]} <directory0/>//<filenums>;<directory1>//<filenums>;... <options>")
ylim = None
titles = ["diff. (order 1)", "Poisson", "int.", "diff. (order 2)", "diff. (order 1 & ana.)", "diff. (order 2 & ana.)"] #, "diff. (order 2)"
types = titles # to plot
for a in sys.argv:
  if a[:6] == "range=": ylim = eval(a[6:])
  elif a[:3] == "hl=": highlighting = eval(a[3:])
  elif a[:6] == "planar":
    planar = True
    direcs += ["../../test_potentials_planar/data/" + folder for folder in planar_direcs]
  elif a[:4] == "dir=":
    direcs = a[4:].split('+')
    if direcs == "planar": direcs = planar_direcs
#  elif a[:7] == "titles=": titles = a[7:].split(',')
#  elif a[:6] == "types=": types = a[6:].split(',')
direc_names = []
potentials = []
area = []
data_left = []
data_right = []
num_add = 3
for directory in direcs:
#  try:
    direc,filenums = directory.split("//")
    filenums = filenums.split(",")
    print(direc)
    #for num in filenums:
    number = 0
    while number < 100:
      num = str(number)
      readout = True
      try:
        maxdev = dat.check_iter(direc, num)
        print(maxdev)
        if maxdev == "nan" or float(maxdev) > 1e-6 or maxdev == "0":
          readout = False
      except: readout = False
      if readout:
        try: pot = dat.get_param("potential", direc+"/")[number]
        except: pot = dat.get_param("potential", direc+"/")
        print(pot, "V")
  #      try:
  #        filename = direc+"/run"+num+"_electrode_charge.dat"
  #        print(filename)
  #        left, right = load_charge(filename)
  #        data_left.append(left+num_add*[-1])
  #        data_right.append(right+num_add*[-1])
  #      except:
        try:
          filename = direc+"/run"+num+"_torder1num_electrode_charge.dat"
          print(filename)
          left, right = load_charge(filename)
          data_left.append(left+num_add*[add])
          data_right.append(right+num_add*[add])
        except:
          break
  #          try:
  #            filename = direc+"/run"+num+"_order1_electrode_charge.dat"
  #            print(filename)
  #            left, right = load_charge(filename)
  #            data_left.append(left+num_add*[-1.])
  #            data_right.append(right+num_add*[-1.])
  #          except:
  #            break
  #      try:
  #        filename = direc+"/run"+num+"_torder2num_electrode_charge.dat"
  #        print(filename)
  #        left, right = load_charge(filename)
  #        data_left[-1][-3] = left[0]
  #        data_right[-1][-3] = right[0]
  #      except:
  #        pass
  #      try:
  #        filename = direc+"/run"+num+"_torder1ana_electrode_charge.dat"
  #        print(filename)
  #        left, right = load_charge(filename)
  #        data_left[-1][-2] = left[0]
  #        data_right[-1][-2] = right[0]
  #      except:
  #        pass
  #      try:
  #        filename = direc+"/run"+num+"_torder2ana_electrode_charge.dat"
  #        print(filename)
  #        left, right = load_charge(filename)
  #        data_left[-1][-1] = left[0]
  #        data_right[-1][-1] = right[0]
  #      except:
  #        pass
  #      try:
  #        filename = direc+"/run"+num+"_order2_electrode_charge.dat"
  #        print(filename)
  #        left, right = load_charge(filename)
  #        data_left[-1][-1] = left[0]
  #        data_right[-1][-1] = right[0]
  #      except:
  #        pass
        potentials.append(pot)
  #      Lx,Ly,Lz = dat.get_param("lengths", direc+"/")
  #      area.append(Lx*Ly)
        direc_names.append(direc+"/")
      number += 1
    print(number)

data_left = np.array(data_left).transpose()
data_right = np.array(data_right).transpose()
#area = np.array(area)

def plotcharge(x, zs):
  for i in range(len(titles)):
    if titles[i] in types:
      if zs[i][0] != add: # -
  #      plt.plot(x, zs[i]/area, dat.markers[i], label=titles[i])
        plt.plot(x, zs[i], dat.markers[i], label=titles[i])
        for h in range(len(x)):
          if direc_names[h] == "corr_test_planar/" or direc_names[h][:5] == "../..":
  #          plt.plot(x[h], zs[i][h]/area[h], 'k'+dat.markers[i])
            plt.plot(x[h], zs[i][h], 'k'+dat.markers[i])
          if highlighting:
            if abs(x[h] - highlight)/highlight < 1e-3:
              gc = dat.get_param("grid_counts", direc_names[h])
              if direc_names[h] == "corr_test_planar/" or direc_names[h][:5] == "../..":
                plt.text(x[h]+1e-3*i, zs[i][h]/area[h], f"planar ({gc[0]}x{gc[1]}x{gc[2]})", color="k") #+direc_names[h])
              else:
                plt.text(x[h]+1e-3*i, zs[i][h]/area[h], f"{gc[0]}x{gc[1]}x{gc[2]}", color=f"c{i}") #+direc_names[h])
    #      for j in range(len(highlight)):
    #        h = highlight[j]
    #        gc = dat.get_param("grid_counts", direcs[h].split("//")[0])
    #        plt.text(x[h], zs[i][h]/area[h], f"{gc[0]}x{gc[1]}x{gc[2]}")
    #        plt.plot(x[h], zs[i][h]/area[h], "k"+dat.markers[i], ms=2)
  if highlighting: plt.legend(ncol=len(titles))
  else: plt.legend()
  if ylim is not None: plt.ylim(*ylim)
  else: plt.ylim(-5,10)
  plt.grid()
  plt.xlabel("potential [V]")
  plt.ylabel("charge [e/nm²]")

if planar == True: suffix = "_all.pdf"
else: suffix = ".pdf"
plotcharge(potentials, data_left)
plt.title("electrode left")
plt.savefig("electrode_charge_left"+suffix)
plt.show()
plt.title("electrode right")
plotcharge(potentials, data_right)
plt.savefig("electrode_charge_right"+suffix)
plt.show()

#direc = "../../potentials/copies/"
#runs = [("run21",[0]),("run25",[1]),("run27",[1])]
#filename = "_3d_profiles.dat"
#titles = []
#plotname = []
#for run in runs:
#  tempT = []
#  tempP = []
#  for num in run[1]:
