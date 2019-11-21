import matplotlib.pyplot as plt  
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from scipy.signal import argrelextrema
     
h = 4.136e-15    # (eV*s)                                           
hbar = 6.582e-16 # (eV*s)
u = 931.49e6     # (eV/c*2) 
m = 10.811*u     # (u)    
c = 2.998e18     # (Angstroem/s)
H = 27.2114      # (eV)

fac_wn = 1./(1.2398e-4)     #  eV to wavenumber

tc = 2.418884*1e-2 # fs = hbar/Ec*1e15
Ec = 27.2114 # 1H = 27.2114 eV

fac_bAA = 0.529177   # 1b = 0.529177 AA

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


R1 = '#e41a1c'
R2 = '#377eb8' 
R3 = '#4daf4a'
R4 = '#984ea3'

# PLOT: MAGNETIZATION ############################################################################## 
file_COOR = open('v0_1.0e-6/td.general/coordinates','r')
COOR = np.loadtxt(file_COOR)
file_COOR.close()

file_E = open('v0_1.0e-6/td.general/energy','r')
E = np.loadtxt(file_E)
file_E.close()

## calculate E_phonon
NN = np.size(COOR[:,1])
Fac = 10

xs_reduced = np.zeros(int(NN/Fac))
data_reduced = np.zeros(int(NN/Fac))

for j in range(int(NN/Fac)):
    xs_reduced[j] =  COOR[j*Fac,1]
    data_reduced[j] = (COOR[j*Fac,3]-COOR[0,3])*fac_bAA
    
POS_MINIMA = argrelextrema(data_reduced, np.greater)[0]
   
WL = 0.0 
T = 0.0  
for ii in range(np.size(POS_MINIMA)-1):
    WL += 1./(np.size(POS_MINIMA)-1)*2*np.pi/(xs_reduced[POS_MINIMA[ii+1]]-xs_reduced[POS_MINIMA[ii]])
    T += 1./(np.size(POS_MINIMA)-1)*(xs_reduced[POS_MINIMA[ii+1]]-xs_reduced[POS_MINIMA[ii]])
 
print(r'T_phon = '+str(T*tc)+' fs') 
print(r'E_phon = '+str(WL*Ec)+' eV')                 # E = hbar*w
print(r'E_phon = '+str(WL*Ec*fac_wn)+' cm^-1')



## PLOT
fig1 = plt.figure(1)
gs1 = gridspec.GridSpec(2, 1)

x_min = 0 
x_max = 100

y_ticks = 3

ax00 = fig1.add_subplot(gs1[0,0]) 
#plt.title(r"$\mathrm{MgB_2}$: $12 \times 12  \times 12$, $\mathrm{Spacing}=0.30$, $\mathrm{E_{2g}}=788.6$ $\mathrm{cm^{-1}}$", y=1.05)
ax00.plot(COOR[:,1]*tc, (COOR[:,3]-COOR[0,3])*fac_bAA, 'k', label=r'$\mathrm{T=}$'+str(np.round(T*tc,2))+'$\mathrm{fs}$')
for m in range(np.size(POS_MINIMA)):
    ax00.axvline(x=xs_reduced[POS_MINIMA[m]]*tc, ymin=-1.0, ymax = +1.0, color='k', linestyle="--", linewidth=1.0)

ax00.set_xlim(x_min, x_max)
ax00.locator_params(nbins=y_ticks, axis='y')
ax00.get_yaxis().set_label_coords(-0.19,0.5)
ax00.set_xticklabels([])
ax00.set_ylabel(r'$\mathrm{u[B_1]}$ ($\mathrm{\AA}$)')


ax10 = fig1.add_subplot(gs1[1,0]) 
ax10.plot(COOR[:,1]*tc, (COOR[:,6]-COOR[0,6])*fac_bAA, 'k')# label=r'$\mathrm{v_0}=1.0 \cdot 10^{-5}$ $\mathrm{bH/\hbar}$')
ax10.set_xlim(x_min, x_max)
ax10.locator_params(nbins=y_ticks, axis='y')
#ax10.get_yaxis().set_label_coords(-0.15,1)
#ax10.set_ylabel('$\mathrm{displacement}$ ($\mathrm{\AA}$)')
#fig1.text(0.045, 0.66, r'$\mathrm{displacement_y}$ $(\mathrm{\AA})$', fontsize=24, ha='center', va='center', rotation='vertical')
ax10.get_yaxis().set_label_coords(-0.19,0.5)
ax10.set_ylabel(r'$\mathrm{u[B_2]}$ ($\mathrm{\AA}$)')
ax10.set_xlabel('$\mathrm{time}$ ($\mathrm{fs}$)')
# =============================================================================
# ax20 = fig1.add_subplot(gs1[2,0]) 
# y_min = +0.9999997
# y_max = +1.0000003
# ax20.plot(E[:,1]*tc, E[:,2], color = R1)
# #ax20.plot(E2[:,1]*tc, E2[:,2]/E2[0,2], color = R2)
# 
# ax20.set_xlabel('$\mathrm{time}$ ($\mathrm{fs}$)')
# ax20.set_ylabel(r'$\mathrm{E_{tot}}/\mathrm{E_{tot}^0}$')
# ax20.set_xlim(x_min, x_max)
# #ax20.set_ylim(y_min, y_max)
# ax20.get_yaxis().set_label_coords(-0.15,0.5)
# ax20.locator_params(nbins=6, axis='x')
# ax20.locator_params(nbins=5, axis='y')
# ax20.get_yaxis().set_label_coords(-0.11,0.5)
# ax20.set_yticks([-9.1, 0.0, 0.1])
# ax20.set_yticklabels([r'$-\mathrm{\pi/2}$', '0' , '$+\mathrm{\pi/2}$'])
# plt.legend(loc='upper left', bbox_to_anchor=(1.00,1.0))
# =============================================================================

plt.subplots_adjust(top=0.95, bottom=0.15, left=0.20, right=0.95, hspace=0.10, wspace=0.05)

#plt.tight_layout()
plt.show()

