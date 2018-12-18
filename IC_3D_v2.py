# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 20:13:08 2017

@author: James Rogers
"""

from __future__ import division
import numpy as np
import numpy.random as rand
import matplotlib.pyplot as plt
import visit_writer
import PM_functions as PM
import Box_class
import Cloud_In_Cell as CIC
import delta_norm


n_s = 0.96               # specral index
h = 0.7                  # cosmological convention
H_0 = (70/3e5)           # Hubble constant in c=1 units
t_0 = 1/H_0              # inverse hubble constant
sigma_8 = 0.9            # sigma_8 cosmological constant today

def Omega_m(a):
    
    omega_m = 1/(a**3)
    
    return omega_m


def transfer_function(k):
    """
    Transfer function for taking modes from radiation
    dominance to matter dominance.
    """
    
    q =  k / (h * h)
    
    x = np.log(1 + 2.34*q) / (2.34*q)
    y = 1 + 3.80*q + (16.2*q)**2 + (5.47*q)**3 + (6.71*q)**4
    
    return x*(y**(-0.25))

def growth_function(a, Omega_M, Omega_Lambda):
    """
    Growth function which describes how modes grow.
    """
    if Omega_M == 1.0 and Omega_Lambda == 0.0:
        return a
    else:
        D_0 = 2.5 * Omega_M * (1/((Omega_M**(4/7)) - Omega_Lambda + (1 + 0.5*Omega_M)*(1 + (1/70)*Omega_M)))    
        return D_0 * a

    
def P(k, As, d3k):
        
    P = np.sqrt(As * (k**(n_s)) * d3k * (transfer_function(k))**2)
    
    return P

def initial_conditions(box, N_g, a, r_0):
    
    L_box = int(N_g*r_0)
    #d3k = (2*np.pi / L_box)**3
    d3k = (24*np.pi**3) / (L_box**3)
    print 'd3k =', d3k
    
    Sx_fft = np.zeros((N_g, N_g, N_g//2 + 1), dtype=np.complex_)
    Sy_fft = np.zeros((N_g, N_g, N_g//2 + 1), dtype=np.complex_)
    Sz_fft = np.zeros((N_g, N_g, N_g//2 + 1), dtype=np.complex_)
    
    Px_fft = np.zeros((N_g, N_g, N_g//2 + 1), dtype=np.complex_)
    Py_fft = np.zeros((N_g, N_g, N_g//2 + 1), dtype=np.complex_)
    Pz_fft = np.zeros((N_g, N_g, N_g//2 + 1), dtype=np.complex_)

    kx = 2 * np.pi * np.fft.fftfreq(N_g, d=r_0) 
    ky = 2 * np.pi * np.fft.fftfreq(N_g, d=r_0) 
    kz = 2 * np.pi * np.fft.rfftfreq(N_g, d=r_0) 

    kmax = 2*np.pi
    D = growth_function(a, box.Omega_M, box.Omega_Lambda)
    print 'D = ',D
    As = D * (sigma_8**2) / delta_norm.sigma_0(kmax)
    print 'As =', As
 
    for x in range(len(kx)):
        for y in range(len(ky)):
            for z in range(len(kz)):
               
                k = np.sqrt((kx[x]*kx[x] + ky[y]*ky[y] + kz[z]*kz[z]))
                
                if k==0:
                    continue
            
                elif k < np.pi:
                    sigma = P(k, As, d3k) / (k**2)
                    a_i = sigma * rand.normal(0,1)
                    b_i = sigma * rand.normal(0,1)
                    c_i = 0.5 * complex(a_i, -b_i)
        
                    Sx_fft[x,y,z] = complex(0., c_i*kx[x])
                    Sy_fft[x,y,z] = complex(0., c_i*ky[y])
                    Sz_fft[x,y,z] = complex(0., c_i*kz[z])
                    
                    Px_fft[x,y,z] = complex(0., c_i*kx[x]*kx[x])
                    Py_fft[x,y,z] = complex(0., c_i*ky[y]*ky[y])
                    Pz_fft[x,y,z] = complex(0., c_i*kz[z]*kz[z])
        
    Sx = (N_g**3) * np.fft.irfftn(Sx_fft) / r_0
    Sy = (N_g**3) * np.fft.irfftn(Sy_fft) / r_0
    Sz = (N_g**3) * np.fft.irfftn(Sz_fft) / r_0
    
    Px = (N_g**3) * np.fft.irfftn(Sx_fft) / r_0
    Py = (N_g**3) * np.fft.irfftn(Sy_fft) / r_0
    Pz = (N_g**3) * np.fft.irfftn(Sz_fft) / r_0
        
    #v_0 = r_0 / t_0
    
    #Px = (N_g**3) * np.fft.irfftn(Px_fft) / v_0
    #Py = (N_g**3) * np.fft.irfftn(Py_fft) / v_0
    #Pz = (N_g**3) * np.fft.irfftn(Pz_fft) / v_0
    
    return Sx,Sy,Sz,Px,Py,Pz


#a_init = 0.01
#da = 0.001
#a = np.arange(a_init,1.001,da)
#N_g = 256
#L_box = 100
#r_0 = L_box/N_g
#
#Omega_M = 1.0
#Omega_Lambda = 0.0
#Omega_K = 0.0
#
#universe = Box_class.Box(N_g, L_box, Omega_M, Omega_K, Omega_Lambda)
#Sx,Sy,Sz,Px,Py,Pz = initial_conditions(universe, N_g, a_init, r_0)
#
##print 'Sx =', Sx
#
#PM.Zeldovich_Approximation(universe,a_init,da,Sx,Sy,Sz,Px,Py,Pz)
#CIC.CIC_density_assignment(universe)
##PM.get_phi(universe, a_init)
#delta = universe.rho - np.ones((N_g,N_g,N_g))
#
#
#X = np.arange(N_g)
#Y = np.arange(N_g)
#x, y = np.meshgrid(X,Y)
#
#plt.figure(1)
#plt.title(r'$\delta$')
#plt.contourf(x, y, delta[8,:,:]) 
#plt.colorbar()
#plt.axis('off')
#plt.show()
#
#plt.figure(2)
#plt.title(r'$S_x$')
#plt.contourf(x, y, Sx[8,:,:]) 
#plt.colorbar()
#plt.axis('off')
#plt.show()
#
#plt.figure(3)
#plt.title(r'$P_x$')
#plt.contourf(x, y, ((a_init-da/2)**2)*Px[8,:,:]/a_init) 
##plt.colorbar()
#plt.axis('off')
#plt.show()

# --------------------------------------------------------------------------- #


#nodal = list(phi_IC.flatten())
#
#zonal = np.ndarray((N_g-1,N_g-1,N_g-1))
#for x in range(N_g-1):
#    for y in range(N_g-1):
#        for z in range(N_g-1):
#            zonal[x,y,z] = (phi_IC[x,y,z] + phi_IC[x+1,y+1,z+1]) / 2
#zonal = list(zonal.flatten())
#"""
#0 for zonal
#1 for nodal
#"""
#
#dims = (N_g, N_g, N_g)
#
#var = (("zonal", 1, 0, zonal), ("nodal", 1, 1, nodal))
#
#visit_writer.WriteRegularMesh("IC.vtk", 0, dims, var)