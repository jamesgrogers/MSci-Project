# -*- coding: utf-8 -*-
"""
Zeldovich approx tests
Created on Mon Oct 30 15:29:49 2017

@author: James Rogers
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations
import datetime
import time
import h5py
import os


import Box_class
import PM_functions as PM
import Cloud_In_Cell as CIC

# ------------------------------------- SETUP ---------------------------------- #

# UNIVERSE SIZE AND MESH
N_g = 16
L_box = 1000
r_0= L_box/N_g

# COSMOLOGICAL PARAMETERS
Omega_M = 1.0
Omega_K = 0.0
Omega_Lambda = 0.0

# INITIALISE
universe = Box_class.Box(N_g, L_box, Omega_M, Omega_K, Omega_Lambda)

a_init = 0.05
a_cross = 10 * a_init
da = 0.005
a = np.arange(a_init, 0.5, da)

current_time_string = datetime.datetime.fromtimestamp(time.time()).strftime('%d.%m.%Y_%H.%M.%S')
newpath = 'ZELDOVICH_OUTPUT/{0}'.format(current_time_string)

if not os.path.exists(newpath):
    os.makedirs(newpath)

A = 1 / ((a_cross) * (2*np.pi/N_g))

print A

PM.zeldovich_test(universe, A, a_init)

#x_history = []
#y_history = []
#z_history = []
#
#px_history = []
#py_history = []
#pz_history = []
#
#rho_history = []
#phi_history = []
#
#x_history.append(list(universe.x))
#y_history.append(list(universe.y))
#z_history.append(list(universe.z))
#px_history.append(list(universe.px))
#py_history.append(list(universe.py))
#pz_history.append(list(universe.pz))    

# ----------------------------------- ITERATE ---------------------------------- #
#PM.VTK_output_x_and_p(universe.x, universe.y, universe.z, universe.px, universe.py, universe.pz, "{0}/universe_{1}.vtk".format(newpath,0))


start_time = time.time()

k = 2*np.pi/N_g

for i in range(1, len(a)):
    
    CIC.CIC_density_assignment(universe)
    #PM.NGP_density_assignment(universe)
    PM.get_phi(universe, a[i])
    CIC.CIC_forward_integrate(universe, a[i], da)  
    #PM.NGP_forward_integrate(universe, a[i], da) 
    
#    x_history.append(list(universe.x))
#    y_history.append(list(universe.y))
#    z_history.append(list(universe.z))
#    px_history.append(list(universe.px))
#    py_history.append(list(universe.py))
#    pz_history.append(list(universe.pz)) 
#    rho_history.append([universe.rho[j,0,0] for j in range(N_g)])
#    phi_history.append([universe.phi[j,0,0] for j in range(N_g)])
    #PM.VTK_output_x_and_p(universe.x,universe.y,universe.z,universe.px,universe.py,universe.pz,"{0}/universe_{1}.vtk".format(newpath,i))


    print 100 * (i)/len(a), '%'
print "--- {0}s seconds ---".format(time.time() - start_time)    
#vx_history = [px_history[i]/a[i] for j in range(len(px_history)-1)]
#vy_history = [py_history[i]/a[i] for j in range(len(py_history)-1)]  
#vz_history = [pz_history[i]/a[i] for j in range(len(pz_history)-1)]
    
#x_history = np.array(x_history)
#y_history = np.array(y_history)
#z_history = np.array(z_history)
#
#px_history = np.array(px_history)
#py_history = np.array(py_history)
#pz_history = np.array(pz_history)
#
#rho_history = np.array(rho_history)
#phi_history = np.array(phi_history)
#
#
#
#
#f = h5py.File("./ZELDOVICH_OUTPUT/{0}/collated_data".format(current_time_string), 'w')
#x_dataset = f.create_dataset('x', data=x_history)
#y_dataset = f.create_dataset('y', data=y_history)
#z_dataset = f.create_dataset('z', data=z_history)
#px_dataset = f.create_dataset('px', data=px_history)
#py_dataset = f.create_dataset('py', data=py_history)
#pz_dataset = f.create_dataset('pz', data=pz_history)
#m_dataset = f.create_dataset('m', data=universe.m)
#a_dataset = f.create_dataset('a', data=a)
#rho_dataset = f.create_dataset('rho', data=rho_history)
#phi_dataset = f.create_dataset('phi', data=phi_history)
#f.close()