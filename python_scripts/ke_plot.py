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
import scipy.integrate as intg

script_path = os.path.dirname(os.path.realpath(__file__))
config_path = pjoin(script_path, '..', 'input.ini')

# Read the configuration file
config = configparser.ConfigParser()
config.read(config_path)

#hdf5 file name and path 
file_name = 'result.h5'
path = sys.argv[1]

path1 = './plots'

path_fig = pjoin(path,path1)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Constants and data loading from input.ini file
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')
NUM_TS = config.getint('time', 'NUM_TS') 
write_interval = config.getint('diagnostics', 'write_interval') 
phase = config.getint('diagnostics', 'write_interval_phase') 
DT_coeff = config.getfloat('diagnostics', 'DT_coeff') 
DATA_TS = int(NUM_TS/write_interval) + 1
t_phase = config.getint('diagnostics', 'write_interval_phase')
DATA_TS_phase = int(NUM_TS /t_phase) + 1

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

n0 = config.getfloat('simulation', 'density')
Te = e*config.getfloat('population', 'tempE')
Ti = e*config.getfloat('population', 'tempI')
Tn = e*config.getfloat('population', 'tempN')
mi = AMU*config.getfloat('population', 'massI')
mn = AMU*config.getfloat('population', 'massN')
nParticlesE = config.getint('population', 'nParticlesE')
EV_TO_K = 11604.52 
#------------------------------------------------------

#------------------------------------------------------
# SIM Vars
NC = config.getint('domain','NC') 
Time = 0
#-----------------------Debye lenght calculation--------
alp = config.getfloat('simulation', 'alpha')
f = config.getfloat('simulation', 'beta')
ni0 = n0
ne0 = n0*((1+alp-f))
ni0 = n0
nn0 = f*ni0  
nb0 = alp*ni0

# Characteristic Debye length
LD = np.sqrt(eps0*Te/(ne0*e**2)) 
electron_spwt = (ne0*NC*LD)/(nParticlesE)
print(electron_spwt)
we = np.sqrt(ne0*e**2/(eps0*me)) # Total Plasma Frequency
#---------------------------------------------------------
#x = np.linspace(0,1025,1025)

#----------------Read hdf5 file ------------
f = h5py.File(pjoin(path, file_name), 'r')

data = f["time_var/kinetic_energy"]
ts = data[:,0]
kee = data[:,1]
kei = data[:,2]
#-------------------------------------------

#-----potential-energy calculation----------
time_steps = sorted(map(int, f['fielddata/efield'].keys()))
PE = np.zeros(len(time_steps))
#print(time_steps)
for i, time_step in enumerate(time_steps):
#for i, time_step in enumerate(ts):
        #time_step = int(time_step/0.05)
        # Get the electric field data for the current time step
        EF_data = f['fielddata/efield/' + str(time_step)]
        #print(len(EF_data))
        x = np.linspace(0,NC,len(EF_data))*LD
        #print((x[1]-x[0])/LD)
        # Integrate the square of the electric field over the spatial domain
        E_sq = (EF_data[:]** 2) * ((me*we**2*LD/e)**2) 
        #print(E_sq)
        integral = intg.trapz(E_sq, x)
        #print(integral)
        # Calculate the potential energy
        PE[i] = 0.5 * eps0 * integral
        #print(PE[i])


THe = ((electron_spwt)*nParticlesE)*Te#*EV_TO_K*kb
#print(THe)
# Normalize the electric potential energy by the total cold electron energy 
PE/= THe
#-------------------------------------------

fig, (ax1,ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(10, 8))


ax1.plot(ts, kee, label='$KE_{e}$')
ax1.set_xlabel('$\omega_{pe}t$')
ax1.set_ylabel('$KE_{e}$')
ax1.grid(True)
ax1.legend(loc='upper right',framealpha=0.5)

ax2.plot(ts, kei, label="$KE_{i}$")
ax2.set_xlabel('$\omega_{pe}t$')
ax2.set_ylabel('$KE_{i}$')
ax2.grid(True)
ax2.legend(loc='upper right',framealpha=0.5)

ax3.plot(ts, PE, label="potential energy")
ax3.set_xlabel('$\omega_{pe}t$')
ax3.set_ylabel('$PE$')
ax3.grid(True)
ax3.legend(loc='upper right',framealpha=0.5)

ax4.plot(ts, kei + kee + PE, label="Total Energy")
ax4.set_xlabel('$\omega_{pe}t$')
ax4.set_ylabel('Total Energy')
ax4.grid(True)
ax4.set_ylim([min(kei + kee + PE) - 1.0, max(kei + kee + PE) + 1.0])
ax4.legend(loc='upper right',framealpha=0.5)

plt.tight_layout()

plt.savefig(pjoin(path,'ke_pe.png'),dpi = 1200)

plt.show()