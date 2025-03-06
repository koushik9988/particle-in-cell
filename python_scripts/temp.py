import numpy as np
import matplotlib.pyplot as plt
from os.path import join as pjoin
import sys
import h5py

# hdf5 file name and path
file_name = 'result.h5'
path = sys.argv[1] #or give location of the path

# Read hdf5 file
f = h5py.File(pjoin(path, file_name), 'r')

metadata_group = f['/metadata']
# Read attributes
NC = metadata_group.attrs['NC']
timestep = 10000
# Plotting
fig, (ax1,ax2,ax3) = plt.subplots(3, 1)

pot = f[f"fielddata/pot/{timestep}"]
efield = f[f"fielddata/efield/{timestep}"]
den_e = f[f"fielddata/den_electron/{timestep}"]
den_i = f[f"fielddata/den_ion/{timestep}"]

x = np.linspace(0, NC, len(pot))
ax1.plot(x, pot[:], color='black', label="Potential")
ax2.plot(x,efield[:], color='blue', label="Electric field")
ax3.plot(x, den_e[:], color='red', label=" electron density")
ax3.plot(x, den_i[:], color='green', label="ion density")
ax1.legend()
ax2.legend()
ax3.legend()
plt.show()
