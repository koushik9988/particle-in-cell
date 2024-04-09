import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import configparser
import h5py


script_path = os.path.dirname(os.path.realpath(__file__))

# Construct the path to the input.ini file
config_path = pjoin(script_path, '..', 'input.ini')

# Read the configuration file
config = configparser.ConfigParser()
config.read(config_path)

file_name = 'result.h5'

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
print(eps0)

EV_TO_K = 11604.52 
tempE = config.getfloat('population', 'tempE')  		     
tempI = config.getfloat('population', 'tempI') 		
tempN = config.getfloat('population', 'tempN') 
tempB = config.getfloat('population', 'tempB') 
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
#-----------------------------------------------------
ni0 = n0
ne0 = n0*((1-alp-f*alp))
nn0 = alp*ni0  
nb0 = f*ni0
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DATA_TS = int(NUM_TS/write_interval) + 1


LD = np.sqrt(eps0*kb*Te*EV_TO_K/(ne0*e**2)) # Characteristic Debye length
we = np.sqrt(ne0*e**2/(eps0*me)) # Total Plasma Frequency


DT = DT_coeff*(1.0/we)

#----------------Read hdf5 file ------------
electric_field_data = []

f = h5py.File(pjoin(path, file_name), 'r')

time_steps = sorted(map(int, f['fielddata/efield'].keys()))

for i, time_step in enumerate(time_steps):
        EF_data = f['fielddata/efield/' + str(time_step)]
        electric_field_data.append(EF_data[:])  # Append electric field data

# Combine electric field data into a 2D array
EF = np.vstack(electric_field_data)

x = np.linspace(0,NC,len(EF.shape[0]))
dx = x[1]-x[0]

print("The shape of EF is: ", EF.shape)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
wpet_1 = 0 #1000
wpet_2 = NUM_TS*DT_coeff
y1 = wpet_1/(DT_coeff*write_interval)
y2 = wpet_2/(DT_coeff*write_interval)

E = EF[int(y1):int(y2),:]
#nde = nde[int(y1):int(y2),:]

#print("The shape of E (reduced EF) is: ", E.shape)
#+++++++++++++++++++++++++ FFT of E-field data ++++++++++++++++++++++++++++++++
F = np.fft.fftn(E, norm='ortho')
#F = np.fft.fftn(nde, norm='ortho')
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define Omega and K
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NUM_TS1 = wpet_1/DT_coeff
NUM_TS2 = wpet_2/DT_coeff
actual_sim_time = (NUM_TS2 - NUM_TS1)*(DT) # Previously it was (NUM_TS*DT)
# -----------------------------------------------------------------------------
omega = 2*np.pi*np.arange(NUM_TS)/(actual_sim_time) #(DT*NUM_TS) #(actual_sim_time) # Unnormalized Omega
k     = 2*np.pi*np.arange(NC+1)/(NC*dx*LD) # Normalized 
print('Length of k: ',len(k))
print('Max of k: ',np.max(k))

Omega, K = np.meshgrid(omega, k, indexing='ij')
print('Shape of Omega: ',Omega.shape)
#------------------------------------------------------------------------------
halflen = np.array(F.shape, dtype=int)//2
Omega = Omega[:halflen[0],:halflen[1]]
K = K[:halflen[0],:halflen[1]]
F = F[:halflen[0],:halflen[1]]
Omega /=we # Normalized Omega
K = K*(LD)  # Normalized K
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++ Raw Analytic Part Calculation ++++++++++++++++++++++++++++++
raw_analytic = True
if raw_analytic:    
    #wla = np.sqrt((we**2*(1 + 0.5*3*k**2*LD**2)))
    wla = np.sqrt(we**2 + (3/2)*(we*LD)**2*k**2)    
    ep = wla/we    # Electron Plasma Wave (Langmuir Wave)
    
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#F_real = np.real(F)
#F_imag = np.imag(F)
Z = np.log(np.abs(F))
#Z = np.abs(F)

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



#fig,ax = plt.subplots(figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
fig, ax = plt.subplots()
#c1 = plt.pcolor(K, Omega, Z,cmap='rainbow',shading='auto',vmin=np.min(Z),vmax=np.max(Z))##
c1 = plt.pcolor(K, Omega, Z, cmap='rainbow',shading='auto',vmin= - 3,vmax = 3)
c1 = plt.contourf(K, Omega, Z, cmap='rainbow',shading='auto',vmin=np.min(Z),vmax=np.max(Z))
#c2 = plt.contourf(K, Omega, z, cmap='rainbow',shading='auto',vmin=np.min(z),vmax=np.max(z))
#plt.contourf(K, Omega, f, cmap='rainbow',shading='auto',vmin=np.min(f),vmax=np.max(f))
#c2 = plt.contour(x*LD,y/we,z , color ='r') ##

print(np.min(Z), np.max(Z))

if raw_analytic:
    plt.plot(k*LD, ep, color='k', linestyle='--', lw = 1.0, label='$EPW$')

# Note: Changing vmin & vmax will change the depth of the plots.
cbar = plt.colorbar()
cbar.set_label('$\zeta$')

ax.set_xlabel('$k \lambda_{D}$')
ax.set_ylabel('$\omega/\omega_{pe}$')
#ax.set_xlim([0, np.max(K)])
ax.set_xlim([0,3.0])
ax.set_ylim([0,5])
leg = ax.legend(loc='upper right',framealpha=0.5)
plt.savefig(pjoin(path,'dispersion.png'),dpi=dpi)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plt.show()

