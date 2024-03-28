"""
Plots the dispersion graph of from the electric field data
Use %reset -f to clear the python workspace
Data File Invoked: processed_results_all.npz
Run as: 
"""
# TO RUN: python3 dispersion_npz.py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import time
import configparser


script_path = os.path.dirname(os.path.realpath(__file__))

# Construct the path to the input.ini file
config_path = pjoin(script_path, '..', 'input.ini')

# Read the configuration file
config = configparser.ConfigParser()
config.read(config_path)
file_name = 'processed_results_all.npz'
#path = './'
path = sys.argv[1]
# ------------------ Comments -------------------------------------------------
# input parameters specific file
# path to the data folder is ../data/data002_vd_20/files for vd=20
# path to the data folder is ../data/data001_vd_80/files for vd=80
#------------------------------------------------------------------------------
# Constants
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
NUM_TS = config.getint('time', 'NUM_TS') 
write_interval = config.getint('diagnostics', 'write_interval') 
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

alp = 0.5
f = 0.9
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
    ndn = data['ndn']
    ndb = data['ndb']
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
ndn = ndn.reshape(DATA_TS, (NC+1))
ndb = ndb.reshape(DATA_TS, (NC+1))
x=x.reshape(DATA_TS,(NC+1))
#x = x[0,:]

#fig, ax = plt.subplots()
#cax = ax.contourf(x, np.arange(DATA_TS), nde, cmap='rainbow',shading='auto')
#fig.colorbar(cax)
#plt.show()
#i=200#DATA_TS-1
#plt.plot(x[i], EF[i],label="E_field")
# Plot each row of the EF matrix
fig, ax = plt.subplots()

"""
for i in range(DATA_TS):
    
    ax.clear()
    i = DATA_TS-1
    ax.plot(x[i], phi[i],label="phi")
    #ax.plot(x[i], EF[i],label="E_field")
    ax.plot(x[i], nde[i],label="nde")
    ax.plot(x[i], ndi[i],label="ndi")
    #ax.plot(x[i], ndn[i],label="ndn")
    #ax.plot(x[i], ndb[i],label="ndb")



    #ax.set_ylim([-30,50])
    #ax.set_ylim([0,5])
    ax.legend()
    plt.pause(1e-8)
"""

i = DATA_TS-1
i = 100
#ax.plot(x[i], phi[i],label="phi")
#ax.plot(x[i], EF[i],label="E_field")
ax.plot(x[i], nde[i],label="nde")
ax.plot(x[i], ndi[i],label="ndi")
ax.plot(x[i], ndn[i],label="ndn")
ax.plot(x[i], ndb[i],label="ndb")
ax.set_ylim([0,2])
ax.legend()
plt.show()