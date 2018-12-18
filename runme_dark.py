"""
Created on Wed Nov 22 17:36:21 2017

@author: James Rogers
"""

from __future__ import division
import numpy as np
import scipy as sp
import datetime
import time
import os
import matplotlib.pyplot as plt
import h5py

import Box_class
import PM_functions as PM
import IC_3D_v2 as IC_3D
import Cloud_In_Cell as CIC

""" --------------------------------- SETUP SAVE LOCATION --------------------------------- """

current_time_string = datetime.datetime.fromtimestamp(time.time()).strftime('%d.%m.%Y_%H.%M.%S')
newpath = 'E:/MSci/CDM/{0}'.format(current_time_string)
if not os.path.exists(newpath):
    os.makedirs(newpath)

""" -------------------- INITIALISE UNIVERSE W/ GIVEN PARAMETERS ------------------------ """

# SCALE FACTOR RANGE
a_init = 0.001
da = 0.001
a = np.arange(a_init,1.001,da)

# UNIVERSE SIZE AND MESH
N_g = 128
L_box = 100
r_0= L_box/N_g

# COSMOLOGICAL PARAMETERS
Omega_M = 1.0
Omega_K = 0.0
Omega_Lambda = 0.0

# INITIALISE
universe = Box_class.Box(N_g, L_box, Omega_M, Omega_K, Omega_Lambda)
PM.save_universe_details(universe,"{0}/Simulation_Details.txt".format(newpath))

"""--------------------------------- READ INITIAL CONDITIONS ------------------------------ """

#IC_file = h5py.File("Initial_Conditions", "r")
#Sx = IC_file["Sx"][:]
#Sy = IC_file["Sy"][:]
#Sz = IC_file["Sz"][:]
#Px = IC_file["Px"][:]
#Py = IC_file["Py"][:]
#Pz = IC_file["Pz"][:]

""" ----------------------------------INITIAL CONDITIONS----------------------------------- """

Sx,Sy,Sz,Px,Py,Pz = IC_3D.initial_conditions(universe, N_g, a_init, r_0)
PM.Zeldovich_Approximation(universe,a_init,da,Sx,Sy,Sz,Px,Py,Pz)
PM.IC_output(Sx,Sy,Sz,Px,Py,Pz,"{0}/Initial_Conditions".format(newpath))
#
""" ------------------------------- OUTPUT FIRST ITERATION -------------------------------- """


#PM.VTK_output_x_and_p(universe.x, universe.y,universe.z,universe.px,universe.py,universe.pz,"{0}/universe_{1}.vtk".format(newpath,0))
PM.HDF5_output_rho(universe.rho,"{0}/rho_{1}".format(newpath,a_init))
PM.HDF5_output_x_and_p(universe.x,universe.y,universe.z,universe.px,universe.py,universe.pz,"{0}/universe_coords_{1}.hdf5".format(newpath,0))

""" ---------------------------- ITERATE PARTICLE MESH METHOD ----------------------------- """
output_counter = 1

for j in range(1, len(a)):

    # CALCULATE DENSITY
    CIC.CIC_density_assignment(universe)

    # CALCULATE GRAV. POTENTIAL
    PM.get_phi(universe, a[j])

    # FORWARD INTEGRATE MOMENTA AND POSITIONS
    CIC.CIC_forward_integrate(universe, a[j], da)

    # I/O
    #if a[j] in [0.2,0.4,0.6,0.8]:

    if a[j] in [0.2,0.4,0.6,0.8]:
        print "RHO OUTPUT"
        PM.HDF5_output_rho(universe.rho,"{0}/rho_{1}".format(newpath,a[j]))
    if j%10 == 0:
        #print 'OUTPUT IN PROGRESS {0}/100'.format(int(j/10))
        #PM.VTK_output_x_and_p(universe.x,universe.y,universe.z,universe.px,universe.py,universe.pz,"{0}/universe_{1}.vtk".format(newpath,j))
        PM.HDF5_output_x_and_p(universe.x,universe.y,universe.z,universe.px,universe.py,universe.pz,"{0}/universe_coords_{1}.hdf5".format(newpath,output_counter))
        output_counter = output_counter + 1
    print 'a =', a[j]

""" --------------------- CALCULATE DENISTY FOR t=t0 POWER SPECTRUM ---------------------- """

CIC.CIC_density_assignment(universe)
PM.HDF5_output_rho(universe.rho,"{0}/rho_1.0".format(newpath))
