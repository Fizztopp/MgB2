import matplotlib.pyplot as plt  
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)
     
h = 4.136e-15    # (eV*s)                                           
hbar = 6.582e-16 # (eV*s)
u = 931.49e6     # (eV/c*2) 
m = 10.811*u     # (u)    
c = 2.998e18     # (Angstroem/s)

tc = 2.418884*1e-2 # fs = hbar/Ec*1e15
Ec = 27.2114 # 1H = 27.2114 eV
fac_wn = 1./(1.2398e-4)     #  eV to wavenumber

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['lines.markersize'] = 5
mpl.rcParams['lines.markeredgewidth'] = 1
mpl.rcParams['font.size'] = 20# <-- change fonsize globally
mpl.rcParams['legend.fontsize'] = 18
mpl.rcParams['axes.titlesize'] = 18
mpl.rcParams['axes.labelsize'] = 18
mpl.rcParams['xtick.major.size'] = 12
mpl.rcParams['ytick.major.size'] = 12
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['xtick.direction'] = 'inout'
mpl.rcParams['ytick.direction'] = 'inout'
mpl.rcParams['figure.titlesize'] = 24
mpl.rcParams['figure.figsize'] = [8.,6.]
mpl.rcParams['text.usetex'] = True

RED = '#e41a1c'
BLUE = '#377eb8'
GREEN = '#4daf4a'
BROWN = '#fdae61'
VIOLETT = '#6a3d9a' 

# PLOT: MAGNETIZATION ############################################################################## 
file_ETOT = open('energies.txt','r')
ETOT = np.loadtxt(file_ETOT)
file_ETOT.close()

file_ETOT_small = open('energies_small.txt','r')
ETOT_small = np.loadtxt(file_ETOT_small)
file_ETOT_small.close()

print(ETOT_small)

a = 3.0391
cc = 3.4866

a1 = a*np.array([0.5,-np.sqrt(3)/2, 0.0])
a2 = a*np.array([0.5,+np.sqrt(3)/2, 0.0])
a3 = cc*np.array([0.0,0.0,1.0])

B1 = 2./3.*a1 + 1./3.*a2 + 0.5*a3 
B2 = 1./3.*a1 + 2./3.*a2 + 0.5*a3 
MG = 0.0

def DELTA(delta):
    B1d = (2./3.-delta)*a1 + (1./3.+delta)*a2 + 0.5*a3
    B10 = 2./3.*a1 + 1./3.0*a2 + 0.5*a3
    return (B1d-B10)[1]

disp_angstroem = np.zeros(np.size(ETOT[:,0]))
for i in range(np.size(ETOT[:,0])):
    disp_angstroem[i] = DELTA(ETOT[i,0])
    print(disp_angstroem[i])

disp_angstroem_small = np.zeros(np.size(ETOT_small[:,0]))
for i in range(np.size(ETOT_small[:,0])):
    disp_angstroem_small[i] = DELTA(ETOT_small[i,0])
    print(disp_angstroem_small[i])

ETOT_small[:,1] = ETOT_small[:,1]*Ec

fit=np.polyfit(disp_angstroem_small, ETOT_small[:,1], 2, rcond=None, full=False, w=None, cov=False)
p = np.poly1d(fit)
pdx2 = np.polyder(p,2)


fig1 = plt.figure(1)
gs1 = gridspec.GridSpec(1, 1)

x_min = np.amin(disp_angstroem)
x_max = np.amax(disp_angstroem)
y_min = np.amin(ETOT[:,1])-np.amin(ETOT[:,1])
y_max = np.amax(ETOT[:,1])-np.amin(ETOT[:,1])

x_plot = np.linspace(np.amin(disp_angstroem_small),np.amax(disp_angstroem_small ),100)

ax00 = fig1.add_subplot(gs1[0,0]) 
#plt.title(r"$\mathrm{MgB_2}$: $12 \times 12  \times 12$, $\mathrm{Spacing}=0.30$, $\mathrm{E_{2g}}=592.6$ $\mathrm{cm^{-1}}$", y=1.02)
ax00.plot(disp_angstroem, ETOT[:,1]-np.amin(ETOT[:,1]), "k-", alpha=0.8, label=r'$\mathrm{DFT}$')
#ax00.plot(ETOT[:,0], pdx2 (ETOT[:,0]), color='GREEN', label=r'FIT')
ax00.set_xlabel('$\mathrm{u[B]}$ $(\mathrm{\AA})$')
ax00.set_ylabel(r'$\mathrm{E_{tot}}$ $(\mathrm{eV})$')
ax00.set_xlim(x_min, x_max)
ax00.set_ylim(y_min, y_max)
#plt.legend()

x_min_s = np.amin(disp_angstroem_small)
x_max_s = np.amax(disp_angstroem_small)
y_min_s = np.amin(ETOT_small[:,1])-np.amin(ETOT_small[:,1])
y_max_s = np.amax(ETOT_small[:,1])-np.amin(ETOT_small[:,1])

ax11 = plt.axes([0,0,1,1])
ip = InsetPosition(ax11, [0.14,0.15,0.8,0.8])

ax11.set_axes_locator(ip)
mark_inset(ax00, ax11, loc1=3, loc2=4, fc="none", ec='0.5')
ax11.plot(disp_angstroem_small, ETOT_small[:,1]-np.amin(ETOT_small[:,1]), "kx", mew=2, label=r'$\mathrm{DFT}$')
#ax11.plot(disp_angstroem, ETOT_small[:,1]-np.amin(ETOT_small[:,1]), "kx", mew=2)
ax11.plot(x_plot, p(x_plot)-np.amin(ETOT_small[:,1]), 'k--', label=r'FIT')
ax11.set_xticks([-0.03, 0.0, 0.03])
ax11.set_xticklabels(['$\mathrm{-0.03}$', '$\mathrm{0}$', '$\mathrm{0.03}$'])
ax00.set_xlim(x_min, x_max)
ax00.set_ylim(y_min, y_max)
ax11.set_xlabel('$\mathrm{u[B]}$ $(\mathrm{\AA})$')
ax11.set_ylabel(r'$\mathrm{E_{tot}}$ $(\mathrm{eV})$')
ax11.locator_params(nbins=4, axis='x')
ax11.locator_params(nbins=4, axis='7')
plt.legend()


MIN=np.argmin(p(x_plot))
print(r'E_phon = '+str(np.sqrt(0.5*pdx2(x_plot[MIN])/(m))*c*hbar)+' eV') ## sqrt(2) from effective mass
print(r'E_phon = '+str(np.sqrt(0.5*pdx2(x_plot[MIN])/(m))*c*hbar*fac_wn)+r' cm^-1') ## sqrt(2) from effective mass


plt.subplots_adjust(top=0.95, bottom=0.15, left=0.12, right=0.95, hspace=0.10, wspace=0.05)

plt.show()

