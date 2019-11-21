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
HH = 27.2114

fac_wn = 1./(1.2398e-4)     #  eV to wavenumber


mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['lines.markersize'] = 5
mpl.rcParams['lines.markeredgewidth'] = 1
mpl.rcParams['font.size'] = 16 # <-- change fonsize globally
mpl.rcParams['legend.fontsize'] = 18
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['xtick.direction'] = 'inout'
mpl.rcParams['ytick.direction'] = 'inout'
mpl.rcParams['figure.titlesize'] = 24
mpl.rcParams['figure.figsize'] = [6.,5]
mpl.rcParams['text.usetex'] = True

RED = '#e41a1c'
BLUE = '#377eb8'
GREEN = '#4daf4a'
BROWN = '#fdae61'
VIOLETT = '#6a3d9a' 

mu = 0.177043
print("$E_F =$ "+str(mu*HH))

file_BANDS = open('static/bandstructure','r')
MAT_BANDS = (np.loadtxt(file_BANDS)-mu)*HH
file_BANDS.close()


N_BAND = np.size(MAT_BANDS[:,0])
k=np.arange(N_BAND)     
    

fig1 = plt.figure(1)
gs1 = gridspec.GridSpec(1, 1)
#fig1.suptitle(r"$\mathrm{MgB_2 (LDA)}$: $12 \times 12  \times 12$, $\mathrm{Spacing}=0.30$", y=1.00)

ax11 = fig1.add_subplot(gs1[0,0])
ax11.set_xticks([0,100, 150, 250, 300, 400])
ax11.set_xticklabels(['$\mathrm{\Gamma}$', '$\mathrm{M}$', '$\mathrm{K}$', '$\mathrm{\Gamma}$', '$\mathrm{A}$', '$\mathrm{L}$'])
ax11.set_ylabel('$\mathrm{Energy}$ $\mathrm{(eV)}$')
ax11.plot(k,MAT_BANDS[:,4], 'k', linewidth=2.0)
ax11.plot(k,MAT_BANDS[:,5], 'k', linewidth=2.0)
ax11.plot(k,MAT_BANDS[:,6], 'k', linewidth=2.0)
ax11.plot(k,MAT_BANDS[:,7], 'k', linewidth=2.0)
ax11.plot(k,MAT_BANDS[:,8], 'k', linewidth=2.0)
ax11.plot(k,MAT_BANDS[:,9], 'k', linewidth=2.0)
ax11.plot(k,MAT_BANDS[:,10], 'k', linewidth=2.0)
ax11.plot(k,MAT_BANDS[:,11], 'k', linewidth=2.0)
ax11.plot(k,MAT_BANDS[:,12], 'k', linewidth=2.0)
ax11.plot(k,MAT_BANDS[:,13], 'k', linewidth=2.0)
ax11.plot(k,MAT_BANDS[:,14], 'k', linewidth=2.0)
ax11.plot(k,[0]*N_BAND, 'k--', linewidth=1.0)
#ax11.plot(k, [mu]*N_BAND, 'k--', label=r"chemical potential")
ax11.set_xlim(0,N_BAND)
ax11.set_ylim(-13,12)
#ax11.grid()
#plt.grid(True)

plt.subplots_adjust(top=0.95, bottom=0.10, left=0.15, right=0.95, hspace=0.15, wspace=0.15)

plt.show()