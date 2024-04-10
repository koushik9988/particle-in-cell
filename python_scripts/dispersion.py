import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import configparser
import h5py


file_name = 'result.h5'
path = sys.argv[1]

#----------------Read hdf5 file ------------
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']


eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

# Read attributes from hdf5 file 
NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval  = metadata_group.attrs['write_int']
write_interval_phase = metadata_group.attrs['write_int_phase']
DT_coeff = metadata_group.attrs['DT_coeff']
nParticlesE = metadata_group.attrs['nE']
nParticlesI = metadata_group.attrs['nI']
nnParticlesN = metadata_group.attrs['nN']
nnParticlesB = metadata_group.attrs['nB']
Te = metadata_group.attrs['Te']
Ti = metadata_group.attrs['Ti']
Tb = metadata_group.attrs['Tb']
Tn = metadata_group.attrs['Tn']
alp = metadata_group.attrs['alpha']
beta = metadata_group.attrs['beta']
mi = metadata_group.attrs['mI']
mn = metadata_group.attrs['mN']
mb = metadata_group.attrs['mB']
n0 = metadata_group.attrs['density']
save_fig = metadata_group.attrs['save_fig']
EV_TO_K = 11604.52 
DATA_TS = int(NUM_TS/write_interval) + 1

#-------------------plasma frequenct and debye lenght calculation---------------------------
ni0 = n0
ne0 = n0*((1-alp-beta*alp))
nn0 = alp*ni0  
nb0 = beta*ni0

LD = np.sqrt(eps0*kb*Te*EV_TO_K/(ne0*e**2)) # Characteristic Debye length
we = np.sqrt(ne0*e**2/(eps0*me)) # Total Plasma Frequency
DT = DT_coeff*(1.0/we)

#----------------------------------------------------------------------------------

#transform elctric field data into a 2D array

electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))

for i, time_step in enumerate(time_steps):
        EF_data = f['fielddata/efield/' + str(time_step)]
        electric_field_data.append(EF_data[:])  # Append electric field data

# Combine electric field data into a 2D array
EF = np.vstack(electric_field_data)

x = np.linspace(0,NC, EF.shape[0])
dx = x[1]-x[0]
print("The shape of EF is: ", EF.shape)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
wpet_1 = 0
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

#-----------plotting --------------------------------------------------------

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

