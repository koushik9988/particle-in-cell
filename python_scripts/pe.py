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
import scipy.integrate as intg
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

file_name = 'processed_results_all.npz'
alp = float(input('Enter alp:'))
file_name_ke = 'ke_%f.txt'%(alp)
#file_name_ke = 'ke_0.800000.txt'
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

EV_TO_K = 11604.52 
tempE = config.getfloat('population', 'tempE')  		     
tempI = config.getfloat('population', 'tempI') 		
tempN = config.getfloat('population', 'tempN') 
tempB = config.getfloat('population', 'tempB') 

nParticlesE = config.getint('population', 'nParticlesE')
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
NUM_TS = config.getint('time', 'NUM_TS') 
write_interval = config.getint('diagnostics', 'write_interval') 
DT_coeff = config.getfloat('diagnostics', 'DT_coeff') 
#------------------------------------------------------
# SIM Vars
NC = config.getint('domain','NC') 
Time = 0
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#vd = int(input('Enter vd:'))
n0 = 1E13
Te = tempE
Ti = tempI
Tn = tempN
Tb = tempB
mi = AMU*config.getint('population', 'massI') 
#mi = 40*AMU
mn = AMU*config.getint('population', 'massN') 
#-------------------------------------------------------
alp = config.getfloat('simulation', 'alpha') 
f = config.getfloat('simulation', 'beta') 
print("value of alpha is :",alp)
print("value of beta is :",f)
#-----------------------------------------------------
ni0 = n0
ne0 = n0*((1-alp-f*alp))
nn0 = alp*ni0  
nb0 = f*ni0
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DATA_TS = int(NUM_TS/write_interval) + 1
LD = np.sqrt(eps0*kb*Te*EV_TO_K/(ne0*e**2)) # Characteristic Debye length
we = np.sqrt(ne0*e**2/(eps0*me)) # Total Plasma Frequency

#LD = np.sqrt(eps0*kb*Ti*EV_TO_K/(ni0*e**2)) # Characteristic Debye length
#we = np.sqrt(ni0*e**2/(eps0*mi)) # Total Plasma Frequency

DT = DT_coeff*(1.0/we)
#++++++++++++++++++++ UNPACK ++++++++++++++++++++++++++++++++++++++++++++++++++
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


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
x = x.reshape(DATA_TS,(NC+1))
x = x[0,:] # SPATIAL GRID: Select only a column, all the columns have the same value
dx = x[1]-x[0] # Normalized StepSize

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
"""
In each row of EF, there are NC+1 = 1025 columns. Data is put in the first row first
and then it goes into the second row. Therefore, all the values of the EF for the
first time step is written in first row in 1025 columns. The data of next time step
is in the second row and so on. Similar is the case for writing and reshaping x.
w_pe*t = NUM_TS*DT_coeff
"""
EF = EF.reshape(DATA_TS,(NC+1))
phi = phi.reshape(DATA_TS,(NC+1))
# create array of the w_pe*t
ts = np.arange(0, NUM_TS+write_interval, write_interval)*DT_coeff

#print("The shape of EF is: ", EF.shape)
xl = NC*dx
electron_spwt = (ne0*xl*LD)/(nParticlesE)
print("electron specific weight ",electron_spwt)
# ++++++++++++++++++++++++++++++  Method-2 : Electric Field ++++++++++++++++++++++++++++++++
# each row of EF correspond to a separate time step. Hence, no of rows implies no of timesteps
PE1 = np.zeros(EF.shape[0])
print(len(PE1))
for i in range(len(PE1)):
    xdata = x*LD # un-normalized x-data   
    ydata = (EF[i,:]**2) * ((me*we**2*LD/e)**2) # un-normalized EF-square data
    # integrate the square electric field over the volume to get the un-normalized electric potential energy
    E_int = intg.trapz(ydata,xdata)        
    # multiply 0.5*eps0 to the un-normalized electric potential energy
    PE1[i] = 0.5*eps0*(E_int)

# Normalize by the thermal energy 
THe = ((electron_spwt)*nParticlesE)*Te*EV_TO_K*kb
# Normalize the electric potential energy by the total cold electron energy 
PE1 /= THe

# +++++++++++++++  Assign Appropriate Potential Energy +++++++++++++++++++++++
Pot_En = PE1  # use Method-2 to calculate Pot_En
#+++++++++++++++++++++ Kinetic Energy ++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++ UNPACK ++++++++++++++++++++++++++++++++++++++++++++++++++
if os.path.exists(pjoin(path,file_name_ke)):
    t, pe, kee, kei, ken, keb = np.loadtxt(pjoin(path,file_name_ke),unpack=True)
else:
    print('No data')
    exit()

ke = kee + kei + ken + keb
#ke /= np.max(ke)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
figsize = np.array([80,80/1.618]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 1200                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24.5 #Screen resolution

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#fig,ax = plt.subplots(3,1, figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
fig,ax = plt.subplots(3,1)
ax[0].plot(ts,Pot_En)
ax[0].set_xlabel('$\omega_{pe}t$')
ax[0].set_ylabel('$PE$')
ax[0].grid(True)

ax[1].plot(t,ke)
ax[1].set_xlabel('$\omega_{pe}t$')
ax[1].set_ylabel('$KE$')
ax[1].grid(True)

ax[2].plot(ts,ke+Pot_En)
ax[2].set_ylim([min(ke+Pot_En)-1.0, max(ke+Pot_En)+1.0])
ax[2].set_xlabel('$\omega_{pe}t$')
ax[2].set_ylabel('$TE$')
ax[2].grid(True)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plt.savefig(pjoin(path,'pe.png'),dpi=dpi)
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plt.show()
