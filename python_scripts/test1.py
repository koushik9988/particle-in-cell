import numpy as np
import matplotlib.pyplot as plt
import os.path
from os.path import join as pjoin
import sys
from scipy.signal import find_peaks

file_name = 'density_fixed.txt'

path = sys.argv[1]

if os.path.exists(pjoin(path, file_name)):
    time, den = np.loadtxt(pjoin(path, file_name), unpack=True)
else:
    print('No data')
    exit()

plt.plot(time, den, label="Electron Density")

# Find peaks in the density data
peaks, _ = find_peaks(den, height=0)

# Plot the peaks
plt.plot(time[peaks], den[peaks], "x", label="Peaks")

# Calculate the time difference between peaks
time_diff = time[peaks[4]] - time[peaks[3]]

# Add labels and legend
plt.xlabel("Time Step")
plt.ylabel("Density")
plt.legend()
plt.grid(True)

# Show the plot
plt.show()

print("Horizontal distance between peaks:", time_diff)
