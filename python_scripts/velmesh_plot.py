import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import time
import configparser
import h5py
import matplotlib.animation as animation

# hdf5 file name and path
file_name = 'result.h5'
path = sys.argv[1]
timestep = int(sys.argv[2]) if len(sys.argv) > 2 else None

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



# Plotting
fig, ax = plt.subplots()

def animate(i):
    j = i * write_interval_phase

    den1 = f[f"fielddata/vel_electron/{j}"]
    den2 = f[f"fielddata/vel_ion/{j}"]
    x = np.linspace(0, NC, len(den1))
    ax.clear()
    ax.plot(x, den1[:], color='black', label="electron vel")
    #ax.plot(x, den2[:], color='red', label="ion vel")
    ax.set_ylabel('$\phi$')
    ax.legend(loc='upper right', framealpha=0.5)

    return ax

if timestep is not None:
    # Show static plot for the specified time step
    j = timestep
    den1 = f[f"fielddata/vel_electron/{j}"]
    den2 = f[f"fielddata/vel_ion/{j}"]
    #den3 = f[f"fielddata/den_negion/{j}"]
    x = np.linspace(0, NC, len(den1))
    ax.plot(x, den1[:], color='black', label="electron vel")
    #ax.plot(x, den2[:], color='red', label="ion vel")
    #ax.plot(x, den3[:], color='green', label="negative ion density")
    ax.legend(loc='upper right', framealpha=0.5)
    plt.show()
else:
    # Animate as usual
    ani = animation.FuncAnimation(fig, animate, frames=DATA_TS_PHASE, blit=False, interval=500, repeat=False)
    plt.show()
