import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import h5py

# hdf5 file name and path
file_name = 'result.h5'
path = sys.argv[1]

path1 = './plots'
path_fig = pjoin(path, path1)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Read hdf5 file
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']

# Constants and data loading from hdf5 file
eps0 = constants('electric constant')
me = constants('electron mass')
e = constants('elementary charge')

# Read individual attributes
NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval  = metadata_group.attrs['write_int']
write_interval_phase = metadata_group.attrs['write_int_phase']
DT_coeff = metadata_group.attrs['DT_coeff']
LD = metadata_group.attrs['LDe']
LDi = metadata_group.attrs['LDi']
wpe = metadata_group.attrs['wpe']
wpi = metadata_group.attrs['wpi']
n0 = metadata_group.attrs['density']
save_fig = metadata_group.attrs['save_fig']
normscheme = metadata_group.attrs['norm_scheme']
numberofspecies = metadata_group.attrs['spno'] 

EV_TO_K = 11604.52

data = f["time_var/kinetic_energy"]

mFactor = wpi/wpe

ts = data[:,0]
PE = data[:,-1]  #last column

if normscheme == 2 or normscheme == 4:
    LD = LDi
    ts *= mFactor

kinetic_data = [data[:, i] for i in range(1, 1 + numberofspecies)] 
kinetic_labels = [f'$KE_{{{i}}}$' for i in range(1, numberofspecies + 1)] 

nrows = 2  
ncols = max(2, numberofspecies)

fig, axes = plt.subplots(nrows, ncols, figsize=(10, 8))

for i in range(numberofspecies):
    ax = axes[0, i] 
    ax.plot(ts, kinetic_data[i], label=kinetic_labels[i])
    ax.set_xlabel('$\omega_{pe}t$')
    ax.set_ylabel(f'{kinetic_labels[i]}')
    ax.grid(True)
    ax.legend(loc='upper right', framealpha=0.5)


ax_pe = axes[1, 0]  
ax_pe.plot(ts, PE, label="Potential Energy")
ax_pe.set_xlabel('$\omega_{pe}t$')
ax_pe.set_ylabel('$PE$')
ax_pe.grid(True)
ax_pe.legend(loc='upper right', framealpha=0.5)


total_energy = np.sum(kinetic_data, axis=0) + PE  
ax_total = axes[1, 1] 
ax_total.plot(ts, total_energy, label="Total Energy")
ax_total.set_xlabel('$\omega_{pe}t$')
ax_total.set_ylabel('Total Energy')
ax_total.set_ylim([min(total_energy) - 1.0, max(total_energy) + 1.0])
ax_total.grid(True)
ax_total.legend(loc='upper right', framealpha=0.5)

plt.tight_layout()

# Save figure 
if save_fig == 1:
    plt.savefig(pjoin(path_fig, 'ke_pe.png'), dpi=1200)

plt.show()
