import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import configparser
import h5py
from scipy.signal import find_peaks, hilbert

# hdf5 file name and path
file_name = 'result.h5'
path = sys.argv[1]

# Open the HDF5 file
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']
NC = metadata_group.attrs['NC']
save_fig = metadata_group.attrs['save_fig']
fielddata_group = f['fielddata/pot']

# Constants
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')
n0 = metadata_group.attrs['density']
write_interval = metadata_group.attrs['write_int']
NUM_TS = metadata_group.attrs['NUM_TS']
ne0 = n0
we = np.sqrt(ne0 * e**2 / (eps0 * me))

# Fixed spatial point
fixed_point_index = 1

data = f["time_var/kinetic_energy"]
ts = data[:, 0]
times = ts * (1 / we)

num_timesteps = int(NUM_TS / write_interval)
print(num_timesteps)
potentials = np.zeros(num_timesteps + 1)

for i in range(num_timesteps + 1):
    j = i * write_interval
    pot = f[f"fielddata/pot/{j}"]
    #pot = f[f"fielddata/efield/{j}"]
    potentials[i] = pot[fixed_point_index]

# Find peaks
peaks, _ = find_peaks(potentials)
peak_times = times[peaks]
peak_diffs = np.diff(peak_times)

fig, ax = plt.subplots()
ax.plot(times, potentials, color='black', label="Potential at fixed point")
ax.set_xlabel("Time")
ax.set_ylabel("Potential")
ax.legend()
plt.savefig(pjoin(path, 'potential_vs_time.png'), dpi=1200)
plt.show()


# Fourier Transform
fft_result = np.fft.fft(potentials)
frequencies = np.fft.fftfreq(len(potentials), d=times[1] - times[0])

# Plot the power spectrum
plt.plot(frequencies, np.abs(fft_result))
plt.title("Frequency Spectrum")
plt.xlabel("Frequency")
plt.ylabel("Amplitude")
plt.savefig(pjoin(path, 'frequency_spectrum.png'), dpi=1200)
plt.show()

# Amplitude Envelope Extraction using Hilbert Transform
analytical_signal = hilbert(potentials)
amplitude_envelope = np.abs(analytical_signal)

plt.plot(times, amplitude_envelope)
plt.title("Amplitude Envelope")
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.savefig(pjoin(path, 'amplitude_envelope.png'), dpi=1200)
plt.show()

# Logarithmic Analysis for Growth/Decay Rate
log_amplitude = np.log(amplitude_envelope)
time = np.arange(len(log_amplitude))
slope, intercept = np.polyfit(time, log_amplitude, 1)

plt.plot(time, log_amplitude, label='Log Amplitude')
#plt.plot(time, slope * time + intercept, label='Fit', linestyle='--')
plt.title("Logarithmic Amplitude Envelope")
plt.xlabel("Time Step")
plt.ylabel("Log Amplitude")
plt.legend()
plt.savefig(pjoin(path, 'log_amplitude_envelope.png'), dpi=1200)
plt.show()

# Print Growth/Decay Rate
print(f"Growth/Decay Rate: {slope}")

# Close the file
f.close()
