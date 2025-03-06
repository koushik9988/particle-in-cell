import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.constants import value as constants
from os.path import join as pjoin
import os
import sys
import matplotlib.animation as animation

# Set path
path = sys.argv[1]#"../data_esw_coll_1e19"

# HDF5 file name and paths for plots
file_name = 'result.h5'
plot_path = './plots'
path_fig = pjoin(path, plot_path)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Read HDF5 file
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']

# Constants
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
LD = metadata_group.attrs['LDe']
LDi = metadata_group.attrs['LDi']
wpe = metadata_group.attrs['wpe']
wpi = metadata_group.attrs['wpi']
save_fig = metadata_group.attrs['save_fig']
normscheme = metadata_group.attrs['norm_scheme']
DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

mfactor = wpi / wpe
if normscheme == 2 or normscheme == 1:
    mfactor = 1
print("mfactor:", mfactor)

data = f["time_var/kinetic_energy"]
ts = data[:, 0]
ts *= mfactor  # converted to w_pi*t

# Get specific weights for electron and beam species
metadata_electron = f.get('/metadata_species/electron')
metadata_beam = f.get('/metadata_species/beam')

if metadata_electron and metadata_beam:
    electron_spwt = metadata_electron.attrs['spwt']
    beam_spwt = metadata_beam.attrs['spwt']
else:
    print("Metadata for 'electron' or 'beam' not found.")
    sys.exit(1)

# Find velocity range across both species
v_min, v_max = float('inf'), float('-inf')
for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase
    data_electron = f[f"particle_electron/{j}"][:, 1]
    data_beam = f[f"particle_beam/{j}"][:, 1]
    combined_data = np.concatenate([data_electron, data_beam])
    v_min = min(v_min, np.min(combined_data))
    v_max = max(v_max, np.max(combined_data))

# Plotting setup
fig, ax = plt.subplots()

# Compute initial histogram (t=0) for combined data with separate scaling
j0 = 0  # Initial time step
data_electron_0 = f[f"particle_electron/{j0}"][:, 1]
data_beam_0 = f[f"particle_beam/{j0}"][:, 1]

nbins = 500
hist_electron_0, bins = np.histogram(data_electron_0, bins=nbins, range=(v_min, v_max), density=True)
hist_beam_0, _ = np.histogram(data_beam_0, bins=nbins, range=(v_min, v_max), density=True)

# Scale each histogram by its specific weight
hist_electron_0 *= electron_spwt
hist_beam_0 *= beam_spwt

# Combine the scaled histograms (average them since they share the same bins)
hist_combined_0 = (hist_electron_0 + hist_beam_0) / 2  # Normalize by averaging
bin_centers_0 = (bins[:-1] + bins[1:]) / 2

def animate(i):
    j = i * write_interval_phase
    data_electron = f[f"particle_electron/{j}"][:, 1]
    data_beam = f[f"particle_beam/{j}"][:, 1]

    # Histogram for species velocities with separate scaling
    nbins = 100
    hist_electron, bins = np.histogram(data_electron, bins=nbins, range=(v_min, v_max), density=True)
    hist_beam, _ = np.histogram(data_beam, bins=nbins, range=(v_min, v_max), density=True)

    # Scale each histogram by its specific weight
    hist_electron *= electron_spwt
    hist_beam *= beam_spwt

    # Combine the scaled histograms
    hist_combined = (hist_electron + hist_beam) / 2  # Average for combined distribution
    bin_centers = (bins[:-1] + bins[1:]) / 2

    # Clear and plot
    ax.clear()

    # Plot the initial combined distribution with dashed line
    ax.plot(bin_centers_0, hist_combined_0, linestyle='--', color='gray', 
            label="Initial Combined (Electron + Beam)")

    # Plot the combined distribution at current time step
    ax.plot(bin_centers, hist_combined, color='blue', 
            label=f"Combined (Electron + Beam) ({path})")
    
    ax.set_xlabel('Velocity (v)')
    ax.set_ylabel('f(v)')
    ax.legend(loc='upper right', framealpha=0.5)
    ax.set_title(f'$\omega_{{pi}}t$ = {ts[i]:.2f} (Combined Distribution)')

    return ax

def on_key(event):
    if event.key == 'enter':
        on_key.frame = min(on_key.frame + 1, DATA_TS_PHASE - 1)
        animate(on_key.frame)
        plt.draw()
    elif event.key == 'backspace':
        on_key.frame = max(on_key.frame - 1, 0)
        animate(on_key.frame)
        plt.draw()
    elif event.key == 'tab':
        on_key.ani = animation.FuncAnimation(fig, animate, frames=DATA_TS_PHASE, blit=False, interval=100, repeat=False)
        plt.draw()

on_key.frame = 0
on_key.ani = None  # To keep reference to animation

fig.canvas.mpl_connect('key_press_event', on_key)

# Start with the first frame
animate(on_key.frame)
plt.show()