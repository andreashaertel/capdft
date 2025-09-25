# Executable used for generating all the plots in ../plots/colloquium

#______________________________________________________________________________
# Imports, settings, definitions
import data_structure as STR
import electrode_charge_density as ECD
import load_data as DAT
import plot_basics as PLOT
import numpy as np
import matplotlib.pyplot as plt
import potentials as POT

# plot settings
plots = "colloquium/" # "../plots/colloquium/"
fontsize = 15
figsize = (3,2)
def marker(i):
  if i==-1: return '+', 'k'
  return '+', f"C{i}"

def draw_walls(left,right):
  plt.axline([left,1],[left,2],color='gray')
  plt.axline([right,1],[right,2],color='gray')
  return

# physical parameters
species = [0,1]
d_range = (0.6,1.4) # range densities

e = 1.602176634e-19 # elementary charge [C]

#______________________________________________________________________________
# data 
# initialize data sets (all potentials for the cdiff-plots and samples at 0.4V)
main = STR.data_set(*STR.get_reduced_dict(STR.main_data, 
            ["potential"]).transpose())
planar = STR.data_set(*STR.get_reduced_dict(STR.planar_data, 
            ["potential"]).transpose())

potentials_m = [0.4]
potentials_p = [0.4]
sample = STR.data_set(*STR.potential_data(potentials_m).transpose(), potentials_m)
planar_sample = STR.data_set(*STR.potential_data(potentials_p, base=STR.planar_data).transpose(), potentials_p)

# load data from files: densities, potentials, electrode charge densities (ECD)
for dat in [sample]:
    dat.get_densities_xz()
    dat.get_total_ES_xz()
for dat in [planar_sample]:
    dat.get_densities_z()
    dat.get_total_ES_z()

methods = ECD.main_ECD_methods
for dat in [main]:
  dat.get_mean_ECD(methods)
methods = ECD.planar_ECD_methods
for dat in [planar]:
  dat.get_mean_ECD(methods)

methods = ECD.main_ECD_methods
for dat in [main]:
  dat.get_local_ECD(methods)

# check that all files have been found
if sample.len != len(potentials_m): print("warning sample")
if planar_sample.len != len(potentials_p): print("warning planar sample")

#______________________________________________________________________________
# Plots
print("Plot:")

# density profiles (structured, heatmaps)
if 1:
  print("density profiles (structured, heatmaps)")
  dat = sample
  for nu,drange in zip(species,[(0.6,1.2),(0.1,1.7)]):
    for d in range(dat.len):
      plt.figure(figsize=(4,1.5))
      plt.axline([0,0],[1,0],c='w',linestyle='--')
      plt.axline([0,1.5],[1,1.5],c='w',linestyle='--')
      if nu==1: # only necessary once
        plt.text(3.2,1.6,'convex',color='w')
        plt.text(3.1,0.1,'concave',color='w')
      POT.heatmap_xz(dat.x_mesh[d],dat.z_mesh[d],dat.densities_xz[nu][d],crange=drange,clabel=r"$\rho$ [nm$^{-3}$]")
      plt.savefig(plots+f"densities_heatmap_nu{nu}.pdf", bbox_inches='tight')
      # plt.show()
d_range

# density profiles (planar)
if 1:
  print("density profiles (planar)")
  for dat in [planar_sample]:
    plt.figure(figsize=figsize)
    draw_walls(0,6.4)
    i = 0
    for nu, label in zip(species,["larger ions", "smaller ions"]):
      for d in range(dat.len):
        z = dat.z(d)
        rho = dat.densities_z(nu,0)[d]
        pot = dat.potentials[d]
        plt.plot(z, rho, '.', color=marker(1-nu)[1], ms=3, label=label)
        i+=1
    PLOT.layout_densities_z()
    plt.legend()
    plt.ylim(*d_range)
    plt.savefig(plots+f"planar_densities.pdf", bbox_inches='tight')
    # plt.show()

# potential profiles (structured)
Lz = 8
if 1:
  print("potential profiles (structured)")
  for pos in [0,0.5]:
    walls = (PLOT.f(pos),Lz-PLOT.f(pos))
    plt.figure(figsize=figsize)
    draw_walls(*walls)
    i = 0
    for dat in [sample]:
      for d in range(dat.len):
        z = dat.z(d)
        vals = dat.total_ES_z(pos)[d]
        pot = dat.potentials[d]
        plt.plot(z, abs(vals), marker(i)[0], color=marker(i)[1], label="DFT")
        if pot > 0 and dat.len <= 3: PLOT.layout_total_ES_z(True, walls, PLOT.volt_reduce_unit(pot), num=i)
        elif pot > 0: PLOT.layout_total_ES_z(False)
        i+=1
    plt.legend(fontsize=8, loc='upper center')
    plt.savefig(plots+f"ES_xpos{pos}.pdf", bbox_inches='tight')
    # plt.show()

# potential profile (planar)
if 1:
  print("potential profile (planar)")
  plt.figure(figsize=figsize)
  draw_walls(0,6.4)
  i = 0
  for dat in [planar_sample]:
    for d in range(dat.len):
      z = dat.z(d)
      vals = np.array(dat.total_ES_z(0)[d])
      pot = dat.potentials[d]
      plt.plot(z, abs(vals), '.', color=marker(i)[1], label="DFT")
      if pot > 0 and dat.len <= 3: PLOT.layout_total_ES_z(True, (0,6.4), PLOT.volt_reduce_unit(pot), num=i)
      elif pot > 0: PLOT.layout_total_ES_z(False)
      i+=1
  plt.legend(fontsize=8)
  plt.savefig(plots+"planar_ES.pdf", bbox_inches='tight')
  # plt.show()

# differential capacitance (planar)
if 1:
  print("differential capacitance (planar)")
  plt.figure(figsize=figsize)
  i = 0
  for side in ["left","right"]:
    ECD.plot_c_diff(planar.potentials, planar.mean_ECD(side), labels=[side+" side"], colors=[f"C{i}."])
    i += 1
  plt.legend(fontsize=8)
  plt.grid()
  plt.savefig(plots+f"mcdiff_planar.pdf", bbox_inches='tight')
  # plt.show()

# differential capacitance (macroscopic mean)
if 1:
  print("differential capacitance (macroscopic mean)")
  plt.figure(figsize=figsize)
  i = 0
  for side in ["left","right"]:
    ECD.plot_c_diff(main.potentials, main.mean_ECD(side), labels=["undulated"], colors=[f"C{i}+"])
    ECD.plot_c_diff(planar.potentials, planar.mean_ECD(side), labels=["planar"], colors=[f"C{i}."])
    i += 1
  plt.legend(ncol=2, title='left:                         right:', fontsize=8, title_fontsize=8, frameon=0)
  plt.grid()
  plt.savefig(plots+f"mcdiff_comp_struct.pdf", bbox_inches='tight')
  # plt.show()

# local differential capacitance
if 1:
  print("local differential capacitance")
  plt.figure(figsize=figsize)
  for s, side in enumerate(["left","right"]):
    ECD.plot_c_diff(planar.potentials, planar.mean_ECD(side), labels=["planar"], 
                    colors=[f"C{s}."], annotate=(side=='left'), annotate_pos = -15)
  dat = main
  methods = ECD.main_ECD_methods
  for pos, region in zip([0,0.5],["concave","convex"]):
    pos = int(pos * len(dat.x_vals[0]))
    ECD_left = [[] for _ in methods]
    ECD_right = [[] for _ in methods]
    for m in range(len(methods)):
      for d in range(dat.len):
        ECD_left[m].append(dat.local_ECD("left")[0][d][pos])
        ECD_right[m].append(dat.local_ECD("right")[0][d][pos])
    ECD.plot_c_diff(dat.potentials,ECD_left,labels=[region],colors=["C0+"], fontsize=PLOT.std_fontsize, annotate=True)
    ECD.plot_c_diff(dat.potentials,ECD_right,labels=[region],colors=["C1+"], fontsize=PLOT.std_fontsize)
  plt.grid()
  plt.ylim(0,13)
  plt.savefig(plots+"local_cdiff.pdf", bbox_inches='tight')
  # plt.show()

#______________________________________________________________________________
# data from Seebeck et al 2022, "Elucidating curvature-capacitance relationships
# in carbon-based supercapacitors" (supplement, fig. S2(a), the graphs labeled
# 'flat' and '\chi_5', extracted via the 'PlotDigitizer')
MD_planar_deltapsi = np.array([-0.002472229432424479, 0.09394308881606209, 0.14585909208894218, 0.19777499351081795, 0.24474663976985342, 0.34610631503218464, 0.4004945477374894, 0.447465990294516, 0.4993818917163916, 0.5488258692588562, 0.5957973118158828, 0.6650185137117167, 0.726823332863291, 0.8207664216793529, 0.85290489504585, 0.8998763376028764, 0.9443756525784829, 0.9938194264189383, 1.0457355315428232, 1.100123560546119, 1.1495673343865742, 1.2113719498361397, 1.248454882067486, 1.3028429110707813, 1.3473424297483967, 1.4017304587516926, 1.4511742325921482, 1.5030901340140237, 1.5599504943007487, 1.6044498092763555, 1.6613101695630796, 1.7033373569572663, 1.7601977172439913, 1.8022249046381775, 1.8788628966822898, 1.9357230532670058, 1.9975278724185799])
MD_planar_cdiff = np.array([2.515100794922934, 2.63590611742561, 2.67114098587421, 2.67114098587421, 2.6157717587577243, 2.4848993606091856, 2.479865719098174, 2.505033719277071, 2.550335663371534, 2.5704698146632605, 2.560402739017398, 2.5100671534119234, 2.454698133671596, 2.4043623406899624, 2.389261623533088, 2.374161113752373, 2.3791947552633843, 2.4043623406899624, 2.434563775003711, 2.4647650019413003, 2.4748322849633224, 2.479865719098174, 2.4848993606091856, 2.5352349462146604, 2.590604173331146, 2.6560404760934953, 2.6761744200090622, 2.6459731930714727, 2.580537097685283, 2.525167663192638, 2.4647650019413003, 2.434563775003711, 2.4194630578468366, 2.4093961895771328, 2.389261623533088, 2.354026962460647, 2.3137584525010353])
MD_chi5_deltapsi = np.array([0, 0.049443773840455715, 0.0988874458299067, 0.1433868626565177, 0.19777499351081795, 0.24969089493269353, 0.296662541191729, 0.34610631503218464, 0.3930779612912201, 0.4499381178759361, 0.4993818917163916, 0.5463535379754271, 0.5982694393973026, 0.6501853408191786, 0.6996291146596342, 0.7515450160815097, 0.7836834894480067, 0.8084055800702435, 0.8504325637624208, 0.8974042100214564, 0.9493203151453413, 0.97651412594498, 1.0531519142870835, 1.100123560546119, 1.1520394619679943, 1.2459825507840567, 1.3028429110707813, 1.3473424297483967, 1.4511742325921482, 1.5006180064326038, 1.5624226218821688, 1.6019776816949354, 1.6514214555353905, 1.7008652293758462, 1.7503090032163018, 1.8071693635030268, 1.8590852649249021, 1.9035845799005084, 1.9555004813223842, 1.9950557448371598])
MD_chi5_cdiff = np.array([2.9278524517245152, 3.058724953561134, 3.038590594893248, 2.8221477426906345, 2.560402739017398, 2.459731567806448, 2.4848993606091856, 2.5302015120798087, 2.5201344364339455, 2.429530340868859, 2.3389262453037727, 2.2785235840524347, 2.2684563010304126, 2.3338928111689206, 2.4647650019413003, 2.6661075517393584, 2.8171143085557824, 2.9278524517245152, 3.104026897655597, 3.179530276063809, 3.073825567029928, 2.947986499328162, 2.580537097685283, 2.3791947552633843, 2.3137584525010353, 2.369127679617521, 2.374161113752373, 2.4093961895771328, 2.6157717587577243, 2.565436380528409, 2.3036913768551726, 2.117449751590148, 1.991610787576461, 2.0167785803791984, 2.1627516956846113, 2.4697988508284707, 2.8573824037630757, 3.1442952002390494, 3.3406041085260965, 3.2449665788261584])

def cdiff_convert_unit(val):
  """from µF/cm² to e/nm²/V"""
  return val * 1e-6/e*1e-14

plt.plot(MD_planar_deltapsi, MD_planar_cdiff, marker(1)[1]+'--', label='planar')
plt.xlabel('$\Delta\psi$ [V]')
plt.ylabel('$c^{diff}$ [µF/cm²]')
plt.plot(MD_chi5_deltapsi, MD_chi5_cdiff, marker(1)[1]+'-',label='undulated')
plt.legend()
# plt.show()

# convert & use additional factor 2 since \delta\psi = 2 \psi^0
print("MD data")
plt.figure(figsize=figsize)
plt.plot(MD_planar_deltapsi/2, 2*cdiff_convert_unit(MD_planar_cdiff), marker(1)[1]+'.', label='planar')
plt.xlabel('$\psi^0$ [V]')
plt.ylabel('$c^{diff}$ [e/nm²/V]')
plt.grid()
plt.savefig(plots+"MD_cdiff_planar.pdf", bbox_inches='tight')
plt.plot(MD_chi5_deltapsi/2, 2*cdiff_convert_unit(MD_chi5_cdiff), marker(1)[1]+'+',label='undulated')
plt.legend(fontsize=8, frameon=0)
plt.savefig(plots+"MD_cdiff_comp_struct.pdf", bbox_inches='tight')
# plt.show()
