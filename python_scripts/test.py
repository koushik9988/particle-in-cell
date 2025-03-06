import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.constants import value as constants
from os.path import join as pjoin
import os
import sys
import matplotlib.animation as animation
import matplotlib as mpp

# Argument
if len(sys.argv) != 3:
    print("Usage: python3 script.py <path> <particle_type>")
    sys.exit(1)

# hdf5 file name and path
file_name = 'result.h5'
path = sys.argv[1]
#particle_type = sys.argv[2]
time = int(sys.argv[2])

plot_path = './plots'

path_fig = pjoin(path, plot_path)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Read hdf5 file
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']

# Constants and data loading from input.ini file
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

# Read attributes
NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval = metadata_group.attrs['write_int']
write_interval_phase = metadata_group.attrs['write_int_phase']
DT_coeff = metadata_group.attrs['DT_coeff']

n0 = metadata_group.attrs['density']
save_fig = metadata_group.attrs['save_fig']

DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1



j = time*write_interval_phase
print(write_interval_phase)


figsize = np.array([100,100/1.618]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 1200                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mpp.rc('text', usetex=False)
mpp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mpp.rc('axes', titlesize=10)
mpp.rc('axes', labelsize=10)
mpp.rc('xtick', labelsize=10)
mpp.rc('ytick', labelsize=10)
mpp.rc('legend', fontsize=10)


# Plotting
fig , ax_phase1 = plt.subplots(figsize=figsize/25.4,constrained_layout=True,dpi=ppi)

data_phase1 = f[f"particle_electron/{j}"]
datax1 = data_phase1[:, 0]
datavx1 = data_phase1[:, 1]

data_phase2 = f[f"particle_ebeam/{j}"]
datax2 = data_phase2[:, 0]
datavx2 = data_phase2[:, 1]


ax_phase1.scatter(datax1, datavx1, marker='.', color='r', alpha= 1, s = 1)
ax_phase1.scatter(datax2, datavx2, marker='.', color='b', alpha= 1, s = 1)
ax_phase1.set_xlabel('$x$')
ax_phase1.set_ylabel('$v$')
ax_phase1.legend(loc='upper right', framealpha=0.5)

plt.savefig(pjoin(path,'phasespace.jpg'),dpi=dpi)   
plt.show()
