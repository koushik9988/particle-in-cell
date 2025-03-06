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


#hdf5 file name and path 
file_name = 'result.h5'
path1 = "../data_esw_coll_1e11"
path2 = "../data_esw_coll_1e15"
path3 = "../data_esw_coll_1e17"
path4 = "../data_esw_coll_1e19"
path5 = "../data_esw_coll_1e20"
path6 = "../data_esw_nocoll"

path7 = './plots'

path_fig = pjoin(path6,path5)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

#----------------Read hdf5 file ------------
f1 = h5py.File(pjoin(path1, file_name), 'r')
f2 = h5py.File(pjoin(path2, file_name), 'r')
f3 = h5py.File(pjoin(path3, file_name), 'r')
f4 = h5py.File(pjoin(path4, file_name), 'r')
f5 = h5py.File(pjoin(path5, file_name), 'r')
f6 = h5py.File(pjoin(path6, file_name), 'r')


data1 = f1["time_var/kinetic_energy"]
data2 = f2["time_var/kinetic_energy"]
data3 = f3["time_var/kinetic_energy"]
data4 = f4["time_var/kinetic_energy"]
data5 = f5["time_var/kinetic_energy"]
data6 = f6["time_var/kinetic_energy"]

ts1 = data1[:,0]
PE1 = data1[:,4]
PE2 = data2[:,4]
PE3 = data3[:,4]
PE4 = data4[:,4]
PE5 = data5[:,4]
PE6 = data6[:,4]


#---------------plotting -------------------------
figsize = np.array([80,80/1.618]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 1200                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)

fig,ax = plt.subplots(figsize= figsize/10.4,constrained_layout=True,dpi=ppi)

ax.plot(ts1, PE1, color ='red',label='elastic collision(1e11)')
ax.plot(ts1, PE2, color ='green',label='elastic collision(1e15)')
ax.plot(ts1, PE3, color ='blue',label='elastic collision(1e17)')
ax.plot(ts1, PE4, color ='purple',label='elastic collision(1e19)')
ax.plot(ts1, PE5, color ='violet',label='elastic collision(1e20)')
ax.plot(ts1, PE6, color = 'black',label='no collision')
ax.set_xlabel('$\omega_{pi}t$')
ax.set_ylabel('Electrostatic potential energy')
#ax.set_ylabel('kbe')
ax.grid(True)
ax.legend(loc='upper right',framealpha=0.5)
#plt.semilogy()
#ax.set_ylim(3e-3,3e1)
plt.tight_layout()

#if(save_fig == 1):
plt.savefig(pjoin(path1,'pe_1e11_1e15_1e17_1e19_1e20_nocoll.png'),dpi = dpi)

plt.show()