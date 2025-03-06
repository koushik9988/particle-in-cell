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




#hdf5 file name and path 
file_name = 'result.h5'

path1 = "../data_run_1"
path2 = "../data_run_2"
path3 = "../data_run_3"
path4 = "../data_run_4"
path5 = "../data_run_5"
path6 = "../data_run_6"
path7 = "../data_run_7"
path8 = "../data_run_8"
path9 = "../data_run_9"
path0 = './plots'


path_fig = pjoin(path1,path0)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

#----------------Read hdf5 file ------------
f1 = h5py.File(pjoin(path1, file_name), 'r')
f2 = h5py.File(pjoin(path2, file_name), 'r')
f3 = h5py.File(pjoin(path3, file_name), 'r')
f4 = h5py.File(pjoin(path4, file_name), 'r')
f5 = h5py.File(pjoin(path5, file_name), 'r')
f6 = h5py.File(pjoin(path6, file_name), 'r')
f7 = h5py.File(pjoin(path7, file_name), 'r')
f8 = h5py.File(pjoin(path8, file_name), 'r')
f9 = h5py.File(pjoin(path9, file_name), 'r')

metadata_group1 = f1['/metadata']
metadata_group2 = f2['/metadata']
metadata_group3 = f3['/metadata']
metadata_group4 = f4['/metadata']
metadata_group5 = f5['/metadata']
metadata_group6 = f6['/metadata']
metadata_group7 = f7['/metadata']
metadata_group8 = f8['/metadata']
metadata_group9 = f9['/metadata']

we1 = metadata_group1.attrs['wpe']
we2 = metadata_group2.attrs['wpe']
we3 = metadata_group3.attrs['wpe']
we4 = metadata_group4.attrs['wpe']
we5 = metadata_group5.attrs['wpe']
we6 = metadata_group6.attrs['wpe']
we7 = metadata_group7.attrs['wpe']
we8 = metadata_group8.attrs['wpe']
we9 = metadata_group9.attrs['wpe']

wp1 = metadata_group1.attrs['wpi']
wp2 = metadata_group2.attrs['wpi']
wp3 = metadata_group3.attrs['wpi']
wp4 = metadata_group4.attrs['wpi']
wp5 = metadata_group5.attrs['wpi']
wp6 = metadata_group6.attrs['wpi']
wp7 = metadata_group7.attrs['wpi']
wp8 = metadata_group8.attrs['wpi']
wp9 = metadata_group9.attrs['wpi']

mfactor1 = wp1/we1
mfactor2 = wp2/we2
mfactor3 = wp3/we3
mfactor4 = wp4/we4
mfactor5 = wp5/we5
mfactor6 = wp6/we6
mfactor7 = wp7/we7
mfactor8 = wp8/we8
mfactor9 = wp9/we9

data1 = f1["time_var/kinetic_energy"]
data2 = f2["time_var/kinetic_energy"]
data3 = f3["time_var/kinetic_energy"]
data4 = f4["time_var/kinetic_energy"]
data5 = f5["time_var/kinetic_energy"]
data6 = f6["time_var/kinetic_energy"]
data7 = f7["time_var/kinetic_energy"]
data8 = f8["time_var/kinetic_energy"]
data9 = f9["time_var/kinetic_energy"]

ts1 = data1[:,0]*mfactor1
ts2 = data2[:,0]*mfactor2
ts3 = data3[:,0]*mfactor3
ts4 = data4[:,0]*mfactor4
ts5 = data5[:,0]*mfactor5
ts6 = data6[:,0]*mfactor6
ts7 = data7[:,0]*mfactor7
ts8 = data8[:,0]*mfactor8
ts9 = data9[:,0]*mfactor9


pe1 = data1[:,5]
pe2 = data2[:,5]
pe3 = data3[:,5]
pe4 = data4[:,5]
pe5 = data5[:,5]
pe6 = data6[:,5]
pe7 = data7[:,5]
pe8 = data8[:,5]
pe9 = data9[:,5]

#ts *= mFactor
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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


fig, ax = plt.subplots(figsize=figsize/10.4,constrained_layout=True,dpi=ppi)

ax.plot(ts1, pe1, label="$\\alpha = 0$", color='cyan') 
ax.plot(ts2, pe2, label="$\\alpha = 0.1$", color='magenta')   
ax.plot(ts3, pe3, label="$\\alpha = 0.5$", color='yellow')
ax.plot(ts4, pe4, label="$\\alpha = 1$", color='black')
ax.plot(ts5, pe5, label="$\\alpha = 3$", color='lime')      
ax.plot(ts6, pe6, label="$\\alpha = 5$", color='orange')
ax.plot(ts7, pe7, label="$\\alpha = 10$", color='purple')
ax.plot(ts8, pe8, label="$\\alpha = 15$", color='red')
ax.plot(ts9, pe9, label="$\\alpha = 20$", color='blue')
ax.set_xlabel("$\\omega_{{pi}} t$")
#ax.set_ylabel(r"$\frac{\frac{1}{2} \epsilon_0 E^2}{Th_{ion}}$")
ax.set_ylabel("potential")
#ax.set_ylabel(r"$\frac{\frac{1}{2} \epsilon_0 E^2}{Th_{ion}}$",fontsize=14)

plt.grid()
ax.legend()
plt.savefig(pjoin(path_fig,'pe_comp_alpha.pdf'),dpi = dpi)
plt.show()
