"""
# This file plots the particle data at fixed times. This is the file to generate plots for the paper.
# Only citing the path to the data file is enough to plot,
# No specific input corresponding to parameteric variation is required
# Run as: 
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import pandas as pd
import seaborn as sns
sns.set(style='whitegrid')
import os.path
from os.path import join as pjoin
import sys
import configparser


script_path = os.path.dirname(os.path.realpath(__file__))

# Construct the path to the input.ini file
config_path = pjoin(script_path, '..', 'input.ini')

# Read the configuration file
config = configparser.ConfigParser()
config.read(config_path)

path = sys.argv[1]
path1 = './files'
path_fig = pjoin(path,path1) #'../PS_figure_new'
# ------------- Comments -------------------------------------------------------
# Path to data file for vd=20 is
# ../data/particle_data/part_vd_20/data
# Path to data file for vd=80 is
# ../data/particle_data/part_vd_80/data
# Figures to be saved in ../data/particle_data/part_vd_20/figs for vd=20 &
# ../data/particle_data/part_vd_80/figs for vd=80
#-------------------------------------------------------------------------------
# Define the time at which data is required
NUM_TS = config.getint('time', 'NUM_TS') 
write_interval = config.getint('diagnostics', 'write_interval') 
DT_coeff = config.getfloat('diagnostics', 'DT_coeff')
v_b = config.getfloat('simulation', 'v_b')
alpha = config.getfloat('simulation', 'alpha')
beta = config.getfloat('simulation', 'beta')

#write_interval = 200
# File index value is k(say), change this index to get the plot at required time.
# Calculate the wpet for a k-value by using DT_coeff beforehand.
k = 50000

wpet = k*DT_coeff # Value of wpet

file_name_ele = 'e%d.txt'%(int(k))
file_name_ion = 'i%d.txt'%(int(k))
file_name_neg = 'n%d.txt'%(int(k))
file_name_beam = 'b%d.txt'%(int(k))

# Load Cold electron data
if os.path.exists(pjoin(path,file_name_ele)):
    xe,ve = np.loadtxt(pjoin(path,file_name_ele),unpack=True, usecols=(0,1))
else:
    print('No data')
    exit()

# Load ion data
if os.path.exists(pjoin(path,file_name_ion)):
    xi,vi = np.loadtxt(pjoin(path,file_name_ion),unpack=True, usecols=(0,1))
else:
    print('No data')
    exit()

# Load negion data
if os.path.exists(pjoin(path,file_name_neg)):
    xn,vn = np.loadtxt(pjoin(path,file_name_neg),unpack=True, usecols=(0,1))
else:
    print('No data')
    exit()

# Load Beam data
if os.path.exists(pjoin(path,file_name_beam)):
    xb,vb = np.loadtxt(pjoin(path,file_name_beam),unpack=True, usecols=(0,1))
else:
    print('No data')
    exit()


# Plot the figure
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
figsize = np.array([150,200]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 1200                          #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24   #Screen resolution

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# In scatter plot s specifies the area
fig,ax = plt.subplots(4,1,figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
# --- Play with the parameter 's' in the plot. For vd=20 at wpet=100, 200 make s = 0.00001 for clarity else use 0.01/0.001

#plt.figure(1)
ax[0].scatter(xe,ve,color='g',marker='.',edgecolor='grey',s=1.0)
ax[0].set_xlabel('$x_{e}$')
ax[0].set_ylabel('$\u03C5_{e}$')
#ax[0].set_ylim(-100, 100)
ax[0].set_title("Electron Phase Space at $\omega_{pe}t = %d$" %(wpet))
#plt.scatter(xe,ve)
#plt.savefig(pjoin(path_fig,'cold_%d_wpet_%d.png'%(v_b, wpet)),dpi=1200)


#plt.figure(2)
ax[1].scatter(xi,vi,color='blue',marker='.',edgecolor='blue',s=1.0)
ax[1].set_xlabel('$x_{i}$')
ax[1].set_ylabel('$\u03C5_{i}$')
#ax[1].set_ylim(-5, 5)
ax[1].set_title("Ion Phase Space at $\omega_{pe}t = %d$" %(wpet))
#plt.savefig(pjoin(path_fig,'cold_%d_wpet_%d.png'%(vb, wpet)),dpi=1200)

#plt.figure(3)
ax[2].scatter(xn,vn,color='r',marker='.',edgecolor='red',s=1.0)
ax[2].set_xlabel('$x_{n}$')
ax[2].set_ylabel('$\u03C5_{n}$')
#ax[2].set_ylim(-0.5,0.5)
ax[2].set_title("Negative Ion Phase Space at $\omega_{pe}t = %d$" %(wpet))
#plt.savefig(pjoin(path_fig,'hot_%d_wpet_%d.png'%(vb, wpet)),dpi=1200)

ax[3].scatter(xb,vb,color='r',marker='.',edgecolor='red',s=1.0)
ax[3].set_xlabel('$x_{b}$')
ax[3].set_ylabel('$\u03C5_{b}$')
#ax[3].set_ylim(-0.5,0.5)
ax[3].set_title("Negative Ion Phase Space at $\omega_{pe}t = %d$" %(wpet))
#plt.savefig(pjoin(path_fig,'hot_%d_wpet_%d.png'%(vb, wpet)),dpi=1200)


# Save Figure
plt.savefig(pjoin(path_fig,'PS_wpet = %d_vb = %f_alpha = %f_beta = %f.png'%(wpet,v_b,alpha,beta)),format='png',dpi=300)

plt.show()
