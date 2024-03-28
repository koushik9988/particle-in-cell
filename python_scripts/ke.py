"""
Plots the kinetic energy of the total system along with cold, hot and beam electrons.
Appropriate files must be chosen using 'file_name' and 'path'.
Use %reset -f to clear the python workspace.
Data File Invoked: ke_1024.txt
Run as:
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import os.path
from os.path import join as pjoin
from scipy.constants import value as constants
import sys
alp = float(input('Enter alp:'))
file_name = 'ke_%f.txt'%(alp)
path = sys.argv[1]
# -------------------------- Comments ------------------------------------------
# input parameters specific file
# path to the data folder is ../data/data002_vd_20/files/ke_1024.txt for vd=20
# path to the data folder is ../data/data001_vd_80/files/ke_1024.txt for vd=80
# vd is an input parameter to run this file and that must be provided along with others.
#++++++++++++++++++++ UNPACK ++++++++++++++++++++++++++++++++++++++++++++++++++
if os.path.exists(pjoin(path,file_name)):
    t,pe,kee,kei,ken,keb = np.loadtxt(pjoin(path,file_name),unpack=True)
else:
    print('No data')
    exit()
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# -------- NORMALIZE THE INDIVIDUAL KINETIC ENERGIES W.R.T THE CODE------------
"""
# In the code velocities were normalized while the mass was not normalized.
# Further, the kinetic energy was converted into eV unit. Hence, here the data
# is deconverted by multiplying by e and divided by me to be a normalized one.
"""
eps0 = constants('electric constant')
e = constants('elementary charge')
me = constants('electron mass')
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
te = kee + kei + ken + keb + pe# Total kinetic energy of the electrons
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
figsize = np.array([200,200/1.618]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 1200                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Plot the figure
fig,ax = plt.subplots(2,2,figsize=figsize/25.4,constrained_layout=True,dpi=ppi)

ax[0][0].plot(t,kee,'r',linestyle='-',linewidth=1.0,label='$KE_{E}$')
#ax[0][0].plot(t,pe,'r',linestyle='-',linewidth=1.0,label='$KE_{E}$')
ax[0][0].set_xlabel('$\omega_{pe}t$')
ax[0][0].set_ylabel('$KE_{E}$')
ax[0][0].legend(loc='upper right',framealpha=0.5)

ax[0][1].plot(t,kei,'g',linestyle='-',linewidth=1.0, label='$KE_{I}$')
ax[0][1].set_xlabel('$\omega_{pe}t$')
ax[0][1].set_ylabel('$KE_{I}$')
ax[0][1].legend(loc='upper right',framealpha=0.5)

ax[1][0].plot(t,ken,'b',linestyle='-',linewidth=1.0, label='$KE_{N}$')
ax[1][0].set_xlabel('$\omega_{pe}t$')
ax[1][0].set_ylabel('$KE_{N}$')
ax[1][0].legend(loc='upper right',framealpha=0.5)

ax[1][1].plot(t,keb,'k',linestyle='-',linewidth=1.0, label='$KE_{B}$')
ax[1][1].set_xlabel('$\omega_{pe}t$')
ax[1][1].set_ylabel('$KE_{B}$')
ax[1][1].legend(loc='upper right',framealpha=0.5)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plt.savefig(pjoin(path,'ke.png'),dpi=dpi)

#fig,ax = plt.subplots(figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
#ax.loglog(t,ke,'b',linestyle='-',linewidth=1.0)

plt.show()
