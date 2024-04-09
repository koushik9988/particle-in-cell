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
write_interval_phase = config.getint('diagnostics', 'write_interval_phase')
DATA_TS_PHASE = int(NUM_TS /write_interval_phase) + 1
save_fig = config.getint('diagnostics', 'save_fig')

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

n0 = config.getfloat('simulation', 'density')
Te = e*config.getfloat('population', 'tempE')
Ti = e*config.getfloat('population', 'tempI')
Tn = e*config.getfloat('population', 'tempN')
mi = AMU*config.getfloat('population', 'massI')
mn = AMU*config.getfloat('population', 'massN')
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
#---------------------------------------------------------
x = np.linspace(0,1025,1025)

#----------------Read hdf5 file ------------
f = h5py.File(pjoin(path, file_name), 'r')
#------------------------------------------

fig, (ax_pot, ax_phase1, ax_phase2, ax_phase3, ax_phase4) = plt.subplots(5, 1, figsize=(10, 8))


def animate(i):
    j = i*write_interval_phase

    data_phase_e = f["particle_electron/%d"%j]
    
    dataex = data_phase_e[:,0]
    dataevx = data_phase_e[:,1]


    data_phase_i = f["particle_ion/%d"%j]
    
    dataix = data_phase_i[:,0]
    dataivx = data_phase_i[:,1]

    data_phase_n = f["particle_negion/%d"%j]
    
    datanx = data_phase_n[:,0]
    datanvx = data_phase_n[:,1]

    data_phase_b = f["particle_beam/%d"%j]
    
    databx = data_phase_b[:,0]
    databvx = data_phase_b[:,1]


    ax_phase1.clear()
    ax_phase1.scatter(dataex,dataevx,marker='.',color='b',alpha=1.0,s=10, label = "electron phase space")
    ax_phase1.set_xlabel('$x$')
    ax_phase1.set_ylabel('$v$')
    ax_phase1.legend(loc='upper right',framealpha=0.5)

    ax_phase2.cla()
    ax_phase2.scatter(dataix,dataivx,marker='.',color='r',alpha=1.0,s=10, label = "ion phase space")
    ax_phase2.set_xlabel('$x$')
    ax_phase2.set_ylabel('$v$')
    ax_phase2.legend(loc='upper right',framealpha=0.5)

    ax_phase3.cla()
    ax_phase3.scatter(datanx,datanvx,marker='.',color='g',alpha=1.0,s=10, label = "negative-ion phase space")
    ax_phase3.set_xlabel('$x$')
    ax_phase3.set_ylabel('$v$')
    ax_phase3.legend(loc='upper right',framealpha=0.5)

    ax_phase4.cla()
    ax_phase4.scatter(databx,databvx,marker='.',color='y',alpha=1.0,s=10, label = "beam phase space")
    ax_phase4.set_xlabel('$x$')
    ax_phase4.set_ylabel('$v$')
    ax_phase4.legend(loc='upper right',framealpha=0.5)

    ax_pot.clear()
    pot = f["fielddata/pot/%d"%j]
    #print(pot[:])
    x = np.linspace(0,NC, len(pot))
    ax_pot.plot(x, pot[:], color='black', label="Potential")
    ax_pot.set_ylabel('$\phi$')
    ax_pot.legend(loc='upper right',framealpha=0.5)

    if(save_fig == 1):
        if(i%1000 == 0):
            plt.savefig(pjoin(path_fig, 'phase-pot_%d.png' % (j)), dpi=1200)

    return ax_phase1, ax_phase2 , ax_phase3, ax_phase4 ,ax_pot


ani = animation.FuncAnimation(fig, animate, frames = DATA_TS_PHASE, blit= False, interval= 100, repeat=False)

plt.show()
