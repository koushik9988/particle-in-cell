import numpy as np 
from scipy.constants import value as constants
import matplotlib as mpp
import matplotlib.pyplot as plt
from sympy import *
# --------------------------------------------------------------------------------------------
AMU = constants('atomic mass constant')
e = constants('elementary charge')
eps = constants('electric constant')
# --------------------------------------------------------------------------------------------
def symbolic_coeffs():
    w, k, vb0, vi0, ve0, vthb, vthi, vthe, vthn, wb, wi, we, wn = symbols('w, k, vb0, vi0, ve0, vthb, vthi, vthe, vthn, wb, wi, we, wn')
    
    A = ((w - k*vb0)**2) - (k**2)*(vthb**2)
    B = ((w - k*vi0)**2) - (k**2)*(vthi**2)
    C = ((w - k*ve0)**2) - (k**2)*(vthe**2)
    D = (k**2)*vthn**2
    expr =  A*B*C*(w**2-D) - (B*C*wb**2)*(w**2-D) - (A*C*wi**2)*(w**2-D) - (A*B*we**2)*(w**2-D) - (A*B*C)*(wn**2)
    p = Poly(expr, w)
    coff = p.all_coeffs()    
    return coff
# --------------------------------------------------------------------------------------------
coff = symbolic_coeffs()
#print(coff)
print('Highest power of omega is: %d\n'%(len(coff)-1))
# --------------------------------------------------------------------------------------------
# Define the densities
alpha = 0.4
beta = 0.1
ni0 = 1E13
ne0 = ni0*(1-alpha- alpha*beta)
nn0 = alpha*ni0
nb0 = beta*nn0
#---------------------------------------------------------------------------------------------
# Define masses 
mi  = 1*AMU
mb = 1*AMU
me = constants('electron mass')
mn  = 1*AMU 
# -------------------------------------------------------------------------------------------
# Define Species Temperatures
Tn = 0.1*e
Tb = 0.1*e
Ti = 0.1*e
Te = 1*e
# Define the polytropic coefficient
gn = 1
gb = 1
gi = 1
ge = 1
# --------------------------------------------------------------------------------------------
# Define the Thermal Velocities
vthb = np.sqrt(gb*Tb/mb) #1E3
vthn = np.sqrt(gn*Tn/mn)
vthi = np.sqrt(gi*Ti/mi)
vthe = np.sqrt(ge*Te/me) #1E4
print(vthi, vthe)
#--------------------------------------------------------------------------------------------
# Define the Drift Velocities
vi0 = 0
vn0 = 0
ve0 = 0
vb0 = -vthe*5
k = np.arange(0, 500, 0.01)
# --------------------------------------------------------------------------------------------
wi = np.sqrt((ni0*e**2)/(eps*mi))
wb = np.sqrt((nb0*e**2)/(eps*mb))
we = np.sqrt((ne0*e**2)/(eps*me))
wn = np.sqrt((nn0*e**2)/(eps*mn))
# --------------------------------------------------------------------------------------------
def EPW_dispersion():
     # Coefficients
     coeff1 = eval(str(coff[0]))
     coeff2 = eval(str(coff[1]))
     coeff3 = eval(str(coff[2]))
     coeff4 = eval(str(coff[3]))
     coeff5 = eval(str(coff[4]))
     coeff6 = eval(str(coff[5]))
     coeff7 = eval(str(coff[6]))
     coeff8 = eval(str(coff[7]))
     coeff9 = eval(str(coff[8]))

     roots = []
     for i in range(1,len(k)): 
         coeffs = [coeff1, coeff2[i], coeff3[i], coeff4[i], coeff5[i], coeff6[i], coeff7[i], coeff8[i], coeff9[i]]
         root = np.roots(coeffs)
         roots.append(root)
     roots = np.array(roots)
     return roots
# --------------------------------------------------------------------------------------------
roots_EPW = EPW_dispersion()
solved_analytic = True
if solved_analytic:    
    ep1 = np.real(roots_EPW[:,0])
    ep2 = np.real(roots_EPW[:,1])
    ep3 = np.real(roots_EPW[:,2])
    ep4 = np.real(roots_EPW[:,3])
    ep5 = np.real(roots_EPW[:,4])
    ep6 = np.real(roots_EPW[:,5])
    ep7 = np.real(roots_EPW[:,6])
    ep8 = np.real(roots_EPW[:,7])

    epim1 = np.imag(roots_EPW[:,0])
    epim2 = np.imag(roots_EPW[:,1])
    epim3 = np.imag(roots_EPW[:,2])
    epim4 = np.imag(roots_EPW[:,3])
    epim5 = np.imag(roots_EPW[:,4])
    epim6 = np.imag(roots_EPW[:,5])
    epim7 = np.imag(roots_EPW[:,6])
    epim8 = np.imag(roots_EPW[:,7])

#print(roots_EPW[:,1])
#print(roots_EPW[:,3])
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
figsize = np.array([80,80/1.618]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 1200                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mpp.rc('text', usetex=False)
mpp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mpp.rc('axes', titlesize=10)
mpp.rc('axes', labelsize=10)
mpp.rc('xtick', labelsize=10)
mpp.rc('ytick', labelsize=10)
mpp.rc('legend', fontsize=10)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if solved_analytic:
    fig, ax = plt.subplots(1, 2, figsize=figsize/10.4,constrained_layout=True,dpi=ppi)
    # ax[0].semilogx(k[1:], ep1/we, color='r', linestyle='-', lw = 2.0, label='$EPW$')
    # ax[0].semilogx(k[1:], ep2/we, color='g', linestyle='-', lw = 2.0, label='$EPW$') 
    # ax[0].semilogx(k[1:], ep3/we, color='b', linestyle='-', lw = 2.0, label='$EPW$') 
    # ax[0].semilogx(k[1:], ep4/we, color='k', linestyle='-', lw = 2.0, label='$EPW$') 
    # ax[0].semilogx(k[1:], ep5/we, color='c', linestyle='-', lw = 2.0, label='$EPW$') 
    # ax[0].semilogx(k[1:], ep6/we, color='m', linestyle='-', lw = 2.0, label='$EPW$')

    ax[0].plot(k[1:], ep1/we, 'r.', markersize = 1.0)
    ax[0].plot(k[1:], ep2/we, 'g.', markersize = 1.0) 
    ax[0].plot(k[1:], ep3/we, 'b.', markersize = 1.0) 
    ax[0].plot(k[1:], ep4/we, 'k.', markersize = 1.0) 
    ax[0].plot(k[1:], ep5/we, 'c.', markersize = 1.0) 
    ax[0].plot(k[1:], ep6/we, 'm.', markersize = 1.0)  
    ax[0].plot(k[1:], ep7/we, 'k.', markersize = 1.0)  
    ax[0].plot(k[1:], ep8/we, 'm.', markersize = 1.0)  
    
    ax[0].set_xlabel('$k \lambda_{D}$')
    ax[0].set_ylabel('$\omega/\omega_{pe}$')   
    ax[0].set_title('Real Roots')
    ax[0].grid(True)
    #ax[0].set_xlim([0, np.max(k)])
    ax[0].set_ylim([0, 2])
    #ax[0].set_xlim([0, 3])

    ax[1].plot(k[1:], epim1/we, 'r.', markersize = 2.0)
    ax[1].plot(k[1:], epim2/we, 'g.', markersize = 2.0)
    ax[1].plot(k[1:], epim3/we, 'b.', markersize = 2.0)
    ax[1].plot(k[1:], epim4/we, 'k.', markersize = 2.0)
    ax[1].plot(k[1:], epim5/we, 'c.', markersize = 2.0)
    ax[1].plot(k[1:], epim6/we, 'm.', markersize = 2.0)
    ax[1].plot(k[1:], epim7/we, 'k.', markersize = 2.0)
    ax[1].plot(k[1:], epim8/we, 'm.', markersize = 2.0)
    
    ax[1].set_xlabel('$k \lambda_{D}$')
    ax[1].set_ylabel('$\omega/\omega_{pe}$')  
    ax[1].set_title('Imaginary Roots') 
    ax[1].grid(True)

plt.show()
