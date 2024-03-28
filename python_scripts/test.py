import numpy as np
import matplotlib.pyplot as plt
import os.path
from os.path import join as pjoin
import sys



file_name = 'density_fixed.txt'

path = sys.argv[1]


if os.path.exists(pjoin(path,file_name)):
    time , den  = np.loadtxt(pjoin(path,file_name),unpack=True)
else:
    print('No data')
    exit()

plt.plot(time, den, label="Electron Density")

# Add labels and legend
plt.xlabel("Time Step")
plt.ylabel("Density")
plt.legend()
plt.grid(True)

# Show the plot
plt.show()
