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
    w, k, v10, v20, ve0, vth1, vth2, vthe, w1, w2, we, wn = symbols('w, k, v10, v20, ve0, vth1, vth2, vthe, w1, w2, we, wn')
    
    A = ((w - k*v10)**2) - (k**2)*(vth1**2)
    B = ((w - k*v20)**2) - (k**2)*(vth2**2)
    C = ((w - k*ve0)**2) - (k**2)*(vthe**2)
    expr =  A*B*C*w**2 - (B*C*w1**2)*(w**2) - (A*C*w2**2)*(w**2) - (A*B*we**2)*(w**2) - (A*B*C)*(wn**2)
    p = Poly(expr, w)
    coff = p.all_coeffs()    
    return coff
# --------------------------------------------------------------------------------------------
coff = symbolic_coeffs()
#print(coff)
print('Highest power of omega is: %d\n'%(len(coff)-1))
# --------------------------------------------------------------------------------------------
# Define the densities
np10 = 0.1E13 # Beam 
np20 = 0.9E13 # background ion
ne0 = 0.5E13
nn0 = 0.5E13
#---------------------------------------------------------------------------------------------
# Define masses 
m1  = 40*AMU
m2 = 2*AMU
me = constants('electron mass')
mn  = 140*AMU 
# -------------------------------------------------------------------------------------------
# Define Species Temperatures
T1 = 0.1*e
T2 = 0.1*e
Te = 1*e
# Define the polytropic coefficient
g1 = 1
g2 = 1
ge = 1
# --------------------------------------------------------------------------------------------
# Define the Thermal Velocities
vth1 = np.sqrt(g1*T1/m1) #1E3
vth2 = np.sqrt(g2*T2/m2)
vthe = np.sqrt(ge*Te/me) #1E4
print(vth1, vthe)
#--------------------------------------------------------------------------------------------
# Define the Drift Velocities
v10 = 1E5
v20 = 0
ve0 = 1E6
k = np.arange(0, 500, 0.01)
# --------------------------------------------------------------------------------------------
w1 = np.sqrt((np10*e**2)/(eps*m1))
w2 = np.sqrt((np20*e**2)/(eps*m2))
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
