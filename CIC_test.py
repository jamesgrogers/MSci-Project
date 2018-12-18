# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 17:43:44 2017

@author: James Rogers
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


import Box_class
import PM_functions as PM
import Cloud_In_Cell as CIC

N_g = 4
#Sx = np.zeros((N_g,N_g,N_g))
#Sy = np.zeros((N_g,N_g,N_g))
#Sz = np.zeros((N_g,N_g,N_g))
#Sx,Sy,Sz = IC_3D.initial_conditions(N_g, 0.01, 1000/N_g)
#print Sx
#print '---'
#print Sy
#print '---'
#print Sz
#print '---'


universe = Box_class.Box(N_g)
universe.m = [1]
universe.x = [0.6]
universe.y = [0.0]
universe.z = [0.0]


CIC.CIC_density_assignment(universe)
print universe.rho


# --------------------------------------------------------------------------- #
#universe = Box_class.Box(N_g)
#Sx,Sy,Sz = IC_3D.initial_conditions(N_g, 0.01, 1000/N_g)
#for i in range(N_g):  
#    
#    for j in range(N_g):  
#        
#        for k in range(N_g):  
#            
#            universe.x.append(i+Sx[i,j,k])
#            universe.y.append(j+Sy[i,j,k])
#            universe.z.append(k+Sz[i,j,k]) 
#universe.m = [1]*N_g*N_g*N_g
#
#CIC.CIC_density_assignment(universe)
#print universe.rho
#plt.figure(3)
#X = np.arange(N_g)
#Y = np.arange(N_g)
#x, y = np.meshgrid(X,Y)
#plt.pcolormesh(x, y, universe.rho[8,:,:]) 
#plt.title('RHO')
#plt.colorbar()
#plt.show()