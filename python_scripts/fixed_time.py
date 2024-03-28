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
path1 = './files/phase_space'

path_fig = pjoin(path,path1) #'../PS_figure_new'

if not os.path.exists(path_fig):
    os.makedirs(path_fig)
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
k = 40000

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
# ... (existing code)

# Plot the figure for Cold Electrons
#fig, ax = plt.subplots(figsize=figsize/25.4, constrained_layout=True, dpi=ppi)
fig, ax = plt.subplots(figsize=figsize/20, dpi=ppi)
ax.scatter(xe, ve, color='g', marker='.', edgecolor='grey', s=1.0)
ax.set_xlabel('$x_{e}$')
ax.set_ylabel('$\u03C5_{e}$')
ax.set_title("Electron Phase Space at $\omega_{pe}t = %d$" % (wpet))
plt.savefig(pjoin(path_fig, 'electron_vb=%d_wpet=%d_alpha=%f_beta=%f.png' % (v_b, wpet,alpha,beta)), dpi=1200)
#plt.show()

# Plot the figure for Ions
fig, ax = plt.subplots(figsize=figsize/25.4, constrained_layout=True, dpi=ppi)
ax.scatter(xi, vi, color='blue', marker='.', edgecolor='blue', s=1.0)
ax.set_xlabel('$x_{i}$')
ax.set_ylabel('$\u03C5_{i}$')
ax.set_title("Ion Phase Space at $\omega_{pe}t = %d$" % (wpet))
plt.savefig(pjoin(path_fig, 'ion_vb=%d_wpet=%d_alpha=%f_beta=%f.png' % (v_b, wpet,alpha,beta)), dpi=1200)
#plt.show()

# Plot the figure for Negative Ions
fig, ax = plt.subplots(figsize=figsize/25.4, constrained_layout=True, dpi=ppi)
ax.scatter(xn, vn, color='r', marker='.', edgecolor='red', s=1.0)
ax.set_xlabel('$x_{n}$')
ax.set_ylabel('$\u03C5_{n}$')
ax.set_title("Negative Ion Phase Space at $\omega_{pe}t = %d$" % (wpet))
plt.savefig(pjoin(path_fig, 'negative_ion_vb=%d_wpet=%d_alpha=%f_beta=%f.png' % (v_b, wpet,alpha,beta)), dpi=1200)
#plt.show()

# Plot the figure for Beams
fig, ax = plt.subplots(figsize=figsize/25.4, constrained_layout=True, dpi=ppi)
ax.scatter(xb, vb, color='r', marker='.', edgecolor='red', s=1.0)
ax.set_xlabel('$x_{b}$')
ax.set_ylabel('$\u03C5_{b}$')
ax.set_title("Beam Phase Space at $\omega_{pe}t = %d$" % (wpet))
plt.savefig(pjoin(path_fig, 'beam_vb=%d_wpet=%d_alpha=%f_beta=%f.png' % (v_b, wpet,alpha,beta)), dpi=1200)
#plt.show()


