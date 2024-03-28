import numpy as np
import matplotlib.pyplot as plt
import os.path
from os.path import join as pjoin
from scipy.constants import value as constants
import sys
import time
import configparser
from matplotlib.animation import FuncAnimation, FFMpegWriter



script_path = os.path.dirname(os.path.realpath(__file__))
config_path = pjoin(script_path, '..', 'input.ini')

# Read the configuration file
config = configparser.ConfigParser()
config.read(config_path)
file_name = 'processed_results_all.npz'
path = sys.argv[1]

# Constants and data loading
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
NUM_TS = config.getint('time', 'NUM_TS') 
write_interval = config.getint('diagnostics', 'write_interval') 
phase = config.getint('diagnostics', 'write_interval_phase') 
DT_coeff = config.getfloat('diagnostics', 'DT_coeff') 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#vd = int(input('Enter vd:'))
n0 = 1E13
Te = 1*e
Ti = 0.1*e
Tn = 0*e
mi = 40*AMU
mn = 140*AMU
#-------------------------------------------------------

#------------------------------------------------------
# SIM Vars
NC = config.getint('domain','NC') 
Time = 0
#-----------------------------------------------------

alp = 0
f = 0
ni0 = n0
ne0 = n0*((1+alp-f))
ni0 = n0
nn0 = f*ni0  
nb0 = alp*ni0
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DATA_TS = int(NUM_TS/write_interval) + 1
print(DATA_TS)

LD = np.sqrt(eps0*Te/(ne0*e**2)) # Characteristic Debye length
#x=np.linspace(0,NC*LD,LD)

if os.path.exists(pjoin(path,file_name)):
    data = np.load(pjoin(path,file_name))
    x = data['x']
    nde = data['nde']
    ndi = data['ndi']
    #ndn = data['ndn']
    #ndb = data['ndb']
    phi = data['phi']
    EF = data['EF']
else:
    print('No data')
    exit()

# Reshape the EF matrix
DATA_TS = int(NUM_TS/write_interval) + 1
EF = EF.reshape(DATA_TS, (NC+1))
phi = phi.reshape(DATA_TS, (NC+1))
nde = nde.reshape(DATA_TS, (NC+1))
ndi = ndi.reshape(DATA_TS, (NC+1))
x=x.reshape(DATA_TS,(NC+1))

# Load phase space plot data
#dir_path = "../data"
dir_path = os.path.dirname(path)

ele = [file for file in os.listdir(dir_path) if file.startswith("e")]
ion = [file for file in os.listdir(dir_path) if file.startswith("i")]

ele = sorted(ele, key=lambda x: int(x[1:-4])) # excluding first letter and txt extension
ion = sorted(ion, key=lambda x: int(x[1:-4]))


# ... (your previous code for phase space plot)

# Plot both dispersion graph and phase space plot on the same graph
t_phase = config.getint('diagnostics', 'write_interval_phase')
DATA_TS = int(NUM_TS / write_interval) + 1
DATA_TS_phase = int(NUM_TS /t_phase) + 1

k = int(DATA_TS/DATA_TS_phase) + 1

#fig, (ax_disp, ax_phase) = plt.subplots(2, 1, figsize=(10, 8))
#fig, (ax_pot, ax_phase1,ax_phase2) = plt.subplots(3, 1, figsize=(10, 8))
fig, (ax_pot, ax_den, ax_phase1, ax_phase2) = plt.subplots(4, 1, figsize=(10, 8))


def update(frame):
    j = frame * k
    ax_pot.clear()
    ax_pot.plot(x[j], phi[j], color='black',label="potential")
    #ax_pot.plot(x[j], EF[j], color='black',label="E_field")
    #ax_pot.plot(x[j], ndi[j], color='red',label="ndi")
    #ax_pot.plot(x[j], nde[j], color='cyan',label="nde")
    ax_pot.set_xlabel("x")
    ax_pot.set_ylabel("$\phi$")
    ax_pot.legend()
    ax_pot.set_title(f"Timestep: {j*write_interval}")

    #density plots
    ax_den.clear()
    ax_den.plot(x[j], ndi[j], color='red',label="ndi")
    ax_den.plot(x[j], nde[j], color='cyan',label="nde")
    ax_den.set_xlabel("x")
    ax_den.set_ylabel("$density$")
    ax_den.legend()
    
    # Plot phase space plot for electron
    ax_phase1.clear()
    with open(os.path.join(dir_path, ele[frame]), 'r') as f:
        data_ele = f.readlines()
    x_ele = [float(line.split()[0]) for line in data_ele]
    v_ele = [float(line.split()[1]) for line in data_ele]
    ax_phase1.scatter(x_ele, v_ele, s=1,color = 'blue', label="electron")
    ax_phase1.set_xlabel("x")
    ax_phase1.set_ylabel("v")
    ax_phase1.legend()
   
    # Set common x-axis for both plots
    ax_phase1.set_xlim(ax_pot.get_xlim())

    # Plot phase space plot for ion
    ax_phase2.clear()
    with open(os.path.join(dir_path, ion[frame]), 'r') as f:
        data_ion = f.readlines()
    x_ion = [float(line.split()[0]) for line in data_ion]
    v_ion = [float(line.split()[1]) for line in data_ion]
    ax_phase2.scatter(x_ion, v_ion, s=1,color= 'red',label="ion")
    #ax_phase2.hist(v_ele, bins= 100,label= "velocity distribution")
    ax_phase2.set_xlabel("x")
    ax_phase2.set_ylabel("v")
    ax_phase2.legend()
   
    # Set common x-axis for both plots
    ax_phase2.set_xlim(ax_pot.get_xlim())

ani = FuncAnimation(fig, update, frames=DATA_TS_phase, interval=50)
# Save the animation as an MP4 file (1080p)
writer = FFMpegWriter(fps= 2, metadata=dict(artist='kaushik'), bitrate=8000)
ani.save('phase_pot.mp4', writer=writer)

