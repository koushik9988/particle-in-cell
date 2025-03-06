import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.constants import value as constants
from os.path import join as pjoin
import os
import sys
import matplotlib.animation as animation

# Argument handling
if len(sys.argv) != 2:
    print("Usage: python3 script.py <path> <particle_type>")
    sys.exit(1)

# Set path and species name from arguments
path = "../data"
path1 = "../data"
species_name = sys.argv[1]

# HDF5 file name and paths for plots
file_name = 'result.h5'
plot_path = './plots'
path_fig = pjoin(path, plot_path)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Read HDF5 file
f = h5py.File(pjoin(path, file_name), 'r')
f1 = h5py.File(pjoin(path1, file_name), 'r')
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
LD = metadata_group.attrs['LDe']
LDi = metadata_group.attrs['LDi']
wpe = metadata_group.attrs['wpe']
wpi = metadata_group.attrs['wpi']
save_fig = metadata_group.attrs['save_fig']
normscheme = metadata_group.attrs['norm_scheme']
DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1


mfactor = wpi/wpe
if normscheme == 2 or normscheme == 1:
    mfactor = 1
print("mfactor:",mfactor)

data = f["time_var/kinetic_energy"]
ts = data[:, 0]
ts *= mfactor  # converted to w_pi*t

# Get specific weight for species if present in metadata
species_metadata = f.get(f'/metadata_species/{species_name}')
if species_metadata:
    species_density = species_metadata.attrs['density']
    species_spwt = 1#(species_density * 1024 * LDi) / 500000
else:
    print(f"Species '{species_name}' not found in metadata.")
    sys.exit(1)

# Find velocity range for the species
v_min, v_max = float('inf'), float('-inf')
for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase
    data_species = f[f"particle_{species_name}/{j}"][:, 1]
    v_min = min(v_min, np.min(data_species))
    v_max = max(v_max, np.max(data_species))

# Plotting setup
fig, ax = plt.subplots()



# Compute initial histogram (t=0)
j0 = 0  # Initial time step
data_species_0 = f[f"particle_{species_name}/{j0}"][:, 1]
data_species1_0 = f1[f"particle_{species_name}/{j0}"][:, 1]

nbins = 500
hist_species_0, bins = np.histogram(data_species_0, bins=nbins, range=(v_min, v_max), density=True)
hist_species1_0, bins = np.histogram(data_species1_0, bins=nbins, range=(v_min, v_max), density=True)

hist_species_0 *= species_spwt  # Scale by specific weight
bin_centers_0 = (bins[:-1] + bins[1:]) / 2

def animate(i):
    j = i * write_interval_phase
    data_species = f[f"particle_{species_name}/{j}"][:, 1]
    data_species1 = f1[f"particle_{species_name}/{j}"][:, 1]

    # Histogram for species velocities
    nbins = 100
    hist_species, bins = np.histogram(data_species, bins=nbins, range=(v_min, v_max), density=True)
    hist_species1, bins = np.histogram(data_species1, bins=nbins, range=(v_min, v_max), density=True)

    # Scale by specific weight
    hist_species *= species_spwt
    bin_centers = (bins[:-1] + bins[1:]) / 2


    

    # Clear and plot
    ax.clear()



    # Plot the initial distribution with dashed/dotted lines
    ax.plot(bin_centers_0, hist_species_0, linestyle='--', color='gray', label=f"Initial {species_name} ({path})")
    #ax.plot(bin_centers_0, hist_species1_0, linestyle='--', color='black', label=f"Initial {species_name} ({path1})")


    ax.plot(bin_centers, hist_species, color='blue', label=f"{species_name.capitalize()} {path}")
    ax.plot(bin_centers, hist_species1, color='red', label=f"{species_name.capitalize()} {path1}")
    #ax.set_title(f'{path}')
    ax.set_xlabel('Velocity (v)')
    ax.set_ylabel('f(v)')
    ax.legend(loc='upper right', framealpha=0.5)
    ax.set_title(f'$\omega_{{pi}}t$ = {ts[i]:.2f} {path}')

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
