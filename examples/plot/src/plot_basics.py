# a few definitions a
import data_structure as STR
import electrode_charge_density as ECD
import load_data as DAT
import numpy as np
import matplotlib.pyplot as plt

std_fontsize = 12

species = [0,1]
e = 1.602176634e-19
kT = 1.380649e-23 * 300
bjerrum = 1.5 #nm
bulk_dens = 0.9 #nm³
debye = 1/np.sqrt(4*np.pi*bjerrum*2*bulk_dens)
print("Debye length: ", debye)
def volt_reduce_unit(voltage):
  return e/kT * voltage

def f(x): return 1.6 * np.sin(np.pi*x)**2

def layout_densities_z(fontsize=std_fontsize):
  plt.grid()
  plt.xlabel("$z$ [nm]", size=fontsize)
  plt.ylabel(r"$\rho$ [nm$^{-3}$]", size=fontsize)
def layout_total_ES_z(PB=False, BCpos=(0,8), BCval=0, num=0, fontsize=std_fontsize):
  plt.grid()
  if PB:
    z = np.linspace(*BCpos)
    phi_left = BCval * np.exp((BCpos[0]-z)/debye)
    phi_right = - BCval * np.exp((z-BCpos[1])/debye)
    plt.plot(z, abs(phi_left), f'C{num}--', label='PB')
    plt.plot(z, abs(phi_right), f'C{num}--')
    plt.ylim(1e-6,1e2)
  plt.semilogy()
  plt.xlabel("$z$ [nm]", size=fontsize)
  plt.ylabel(r"$|\phi|$ [kT]", size=fontsize)
def layout_local_ECD(fontsize=std_fontsize):
  plt.grid()
  plt.xlabel("$x$ [nm]", size=fontsize)
  plt.ylabel(r"$\sigma$ [e/nm²]", size=fontsize)
def layout_mean_ECD(fontsize=std_fontsize):
  plt.grid()
  plt.xlabel("$\psi$ [V]", size=fontsize)
  plt.ylabel(r"$\hat\sigma$ [e/nm²]", size=fontsize)

