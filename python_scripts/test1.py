import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import configparser
import h5py
from sympy import *


file_name = 'result.h5'
path = sys.argv[1]

#----------------Read hdf5 file ------------
#f = h5py.File(pjoin(path, file_name), 'r', driver='core', map='read')
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']
metadata_electron = f['/metadata_species/electron']
metadata_ion = f['/metadata_species/ion']


eps0 = constants('electric constant')
eps = constants('electric constant')
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
LD = metadata_group.attrs['LDe']
LDi = metadata_group.attrs['LDi']
we = metadata_group.attrs['wpe']
wp = metadata_group.attrs['wpi']
n0 = metadata_group.attrs['density']
save_fig = metadata_group.attrs['save_fig']
normscheme = metadata_group.attrs['norm_scheme']

Te = metadata_electron.attrs['temperature']
Ti = metadata_ion.attrs['temperature']
mi = metadata_ion.attrs['mass']

EV_TO_K = 11604.52 
DATA_TS = int(NUM_TS/write_interval) + 1

DT = DT_coeff*(1.0/we)
metadata_negion = f['/metadata_species/negion']
metadata_beam = f['/metadata_species/beam']

alp = metadata_negion.attrs["density"]
beta = metadata_beam.attrs["density"]
#----------------------------------------------------------------------------------
#transform elctric field data into a 2D array

electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))

for i, time_step in enumerate(time_steps):
        EF_data = f['fielddata/efield/' + str(time_step)]
        electric_field_data.append(EF_data[:])  # Append electric field data

# Combine electric field data into a 2D array
EF = np.vstack(electric_field_data)

x = np.linspace(0,NC, EF.shape[1])
dx = x[1]-x[0]
print("The shape of EF is: ", EF.shape[1])
print("The value of dx is: ", dx)

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
Omega /= we # Normalized Omega
K = K*(LD)  # Normalized K
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++ Raw Analytic Part Calculation ++++++++++++++++++++++++++++++
raw_analytic_EPW = False
if raw_analytic_EPW:    
    #epe = np.sqrt((we**2*(1 + 0.5*3*k**2*LD**2)))
    epe = np.sqrt(we**2 + 2*(3/2)*(we*LD)**2*k**2)
    #epe = np.sqrt(we**2)# + 2*(3/2)*(we*LD)**2*k**2)       
    #epe = wla/we    # Electron Plasma Wave (Langmuir Wave)


raw_analytic_PW = False
if raw_analytic_PW:    
    pw = np.full(len(k), we)
    #epe = epe/we


raw_analytic_IAW = False
if raw_analytic_IAW:    
    epi = np.sqrt((Te*EV_TO_K*kb/mi)*(1/(1+k**2*LD**2)) + (Ti*EV_TO_K*kb/mi))*k 
    #epi = np.sqrt((Te*EV_TO_K*kb/mi)*(1/(1+k**2*LDi**2)))*k
    #epi = wla/wp    # Electron Plasma Wave (Langmuir Wave)

alpha1 = alp/(1+alp+beta) #0#alp/(1+alp+beta)
raw_analytic_ILW = True
if raw_analytic_ILW:    
    ilw = np.sqrt(((Ti*EV_TO_K*kb/mi)*(1 + alpha1 + 3*(k**2*LDi**2) + 3*(1-alpha1)*(Ti/Te)))/((k*LDi)**2 + (Ti/Te)*(1-alpha1)))*k


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def symbolic_coeffs():
    #w, k, v10, v20, ve0, vth1, vth2, vthe, w1, w2, we, wn = symbols('w, k, v10, v20, ve0, vth1, vth2, vthe, w1, w2, we, wn')
    w, k, vb ,vthe, we, wi, wn, wb = symbols('w, k, vb ,vthe, we, wi, wn, wb')
    
    A = k**2*vthe**2
    B = w**2
    C = ((w - k*vb)**2)
    expr =  A*B*C + we**2*B*C - wi**2*A*C - wn**2*A*C - wb**2*A*B
    p = Poly(expr, w)
    coff = p.all_coeffs()    
    return coff
# --------------------------------------------------------------------------------------------
coff = symbolic_coeffs()
#print(coff)
print('Highest power of omega is: %d\n'%(len(coff)-1))

ni0 = metadata_group.attrs['density']
alpha = alp
beta = beta
ne0 = ni0/(1+alpha+beta)
nn0 = alpha*ne0
nb0 = beta*ne0

# Define masses 
me = constants('electron mass')
mi = 1*AMU
mb = 1*AMU
mn = 1*AMU 

vthe = np.sqrt(Te*e/me) #1E4
vthi = np.sqrt(Ti*e/mi)
vthb = vthi
vthn = vthi

vb = 10*vthi

we = np.sqrt((ne0*e**2)/(eps*me))
wn = np.sqrt((nn0*e**2)/(eps*mn))
wi = np.sqrt((ni0*e**2)/(eps*mn))
wb = np.sqrt((nb0*e**2)/(eps*mn))

def EPW_dispersion():
     # Coefficients
     coeff1 = eval(str(coff[0]))
     coeff2 = eval(str(coff[1]))
     coeff3 = eval(str(coff[2]))
     coeff4 = eval(str(coff[3]))
     coeff5 = eval(str(coff[4]))
   

     roots = []
     for i in range(1,len(k)): 
         coeffs = [coeff1[i], coeff2[i], coeff3[i], coeff4[i], coeff5[i]]
         root = np.roots(coeffs)
         roots.append(root)
     roots = np.array(roots)
     return roots
# -------------------

roots_EPW = EPW_dispersion()

solved_analytic = True
#LD = LDi

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



ep1 = np.real(roots_EPW[:,0])
ep2 = np.real(roots_EPW[:,1])
ep3 = np.real(roots_EPW[:,2])
ep4 = np.real(roots_EPW[:,3])
#ep5 = np.real(roots_EPW[:,4])


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#F_real = np.real(F)
#F_imag = np.imag(F)
Z = np.log(np.abs(F))
#Z = np.abs(F)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#-----------plotting --------------------------------------------------------

#fig,ax = plt.subplots(figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
fig, ax = plt.subplots()
#c1 = plt.pcolor(K, Omega, Z,cmap='rainbow',shading='auto',vmin=np.min(Z),vmax=np.max(Z))##
c1 = plt.pcolor(K, Omega, Z, cmap='rainbow',shading='auto',vmin= -7,vmax = 7)
#c1 = plt.contourf(K, Omega, Z, cmap='rainbow',shading='auto',vmin=np.min(Z),vmax=np.max(Z))

print(np.min(Z), np.max(Z))


if normscheme == 2:
    LD = LDi
    we = wp
if normscheme == 4:
    LD = LDi
    we = we    
if normscheme == 1 or normscheme == 3:
    LD = LD
    we = we


if raw_analytic_EPW:
    plt.plot(k*LD, epe/we, color='k', linestyle='--', lw = 1.0, label='$EPW$')
    #plt.plot(k*LD, ep, color='k', linestyle='--', lw = 1.0, label='$IAW$')

if raw_analytic_PW:
     plt.plot(k*LD, pw/we, color='k', linestyle='--', lw = 1.0, label='$plasma oscillation$')


if raw_analytic_IAW:
    #plt.plot(k*LD, ep, color='k', linestyle='--', lw = 1.0, label='$EPW$')
    plt.plot(k*LD, epi/we, color='k', linestyle='dotted', lw = 1.0, label='$IAW$')

if raw_analytic_ILW:
    #plt.plot(k*LD, ep, color='k', linestyle='--', lw = 1.0, label='$EPW$')
    plt.plot(k*LD, ilw/we, color='k', linestyle='-.', lw = 1.0, label='$IAW$')
    plt.plot(k[1:]*LD, ep1/we, 'r.', markersize = 1.0)
    plt.plot(k[1:]*LD, ep2/we, 'g.', markersize = 1.0) 
    plt.plot(k[1:]*LD, ep3/we, 'b.', markersize = 1.0) 
    plt.plot(k[1:]*LD, ep4/we, 'k.', markersize = 1.0) 
    #plt.plot(k[1:]*LD, ep5/we, 'c.', markersize = 1.0) 

# Note: Changing vmin & vmax will change the depth of the plots.e
cbar = plt.colorbar()
cbar.set_label('$\zeta$')

ax.set_xlabel('$k \lambda_{D}$')
ax.set_ylabel('$\omega/\omega_{pe}$')
#ax.set_xlim([0, np.max(K)])
ax.set_xlim([0,3])
ax.set_ylim([0,1.5])
leg = ax.legend(loc='upper right',framealpha=0.5)
plt.savefig(pjoin(path,'dispersion.pdf'),dpi=dpi)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plt.show()

