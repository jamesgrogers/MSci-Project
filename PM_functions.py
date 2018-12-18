"""
MSci Project - PM code
@author: James Rogers

All dimensionless variables follow that of Anatoly Klypin's PM code 
"""
from __future__ import division
import os
import numpy as np
from numpy import fft
import h5py
import sys
sys.path.insert(0, 'C:\Users\James Rogers\AppData\Local\Programs\LLNL\VisIt 2.13.0\lib')
import visit_writer

Omega_0 = 1.0

# ------------------------------------ INITIAL PARTICLE POSITION ASSIGNMENT ----------------------------- #

def Zeldovich_Approximation(box, a_init, da, Sx, Sy, Sz, Px, Py, Pz):
    
    # perturb each particle from cell centre using displacement arrays from Zel'dovich step
    # also perturb each particle momentum using similar arrays (See IC_3D.py)
    
    for i in range(0,box.N_g):  
    
        for j in range(0,box.N_g):  
        
            for k in range(0,box.N_g):  
            
                box.x.append(i+Sx[i,j,k])
                box.y.append(j+Sy[i,j,k])
                box.z.append(k+Sz[i,j,k]) 
                
                box.px.append(-((a_init-(da/2))**2)*Px[i,j,k]/a_init)
                box.py.append(-((a_init-(da/2))**2)*Py[i,j,k]/a_init)
                box.pz.append(-((a_init-(da/2))**2)*Pz[i,j,k]/a_init)
                 
    # Periodic Boundary conditions            
                
    for i in range(len(box.x)):
        if box.x[i] >= box.N_g:
            box.x[i] = box.x[i] % box.N_g
        if box.x[i] < 0:
            box.x[i] = box.N_g - (np.abs(box.x[i]) % box.N_g)
    
        if box.y[i] >= box.N_g:
            box.y[i] = box.y[i] % box.N_g
        if box.y[i] < 0:
            box.y[i] = box.N_g - (np.abs(box.y[i]) % box.N_g)
    
        if box.z[i] >= box.N_g:
            box.z[i] = box.z[i] % box.N_g
        if box.z[i] < 0:
            box.z[i] = box.N_g - (np.abs(box.z[i]) % box.N_g)
    return



def zeldovich_test(box, A, a_init):
    
    k = 2*np.pi / box.N_g
    
    for i in range(0,box.N_g):  # y direction
    
        for j in range(0,box.N_g):  # z direction
    
            for q in range(0,box.N_g):  # x direction
            
                box.z.append(j)
                box.y.append(i)
                box.x.append(q + a_init*A*np.sin(k*q))
                
                box.px.append(a_init*a_init*A*np.sin(k*q))
                box.py.append(0.0)
                box.pz.append(0.0)
    

            

# --------------------------------------------- DENSITY ASSIGNMENT -------------------------------------- #

def NGP_density_assignment(box):
    """
    Nearest grid point density assignment assumes particles are point like and all of particle's mass is 
    assigned to single grid cell that contains it
    """    
    
    box.rho = np.zeros((box.N_g, box.N_g, box.N_g))
    
    for i in range(len(box.x)):
        x_cell = int(box.x[i] // 1)
        y_cell = int(box.y[i] // 1)
        z_cell = int(box.z[i] // 1)
        
        
        if x_cell >= box.N_g:
            x_cell = x_cell % box.N_g
        
        if y_cell >= box.N_g:
            y_cell = y_cell % box.N_g
        
        if z_cell >= box.N_g:
            z_cell = z_cell % box.N_g
        
        box.rho[x_cell,y_cell,z_cell] += box.m[i]
    
    # # # # # ONLY NECESSARY FOR ZELDOVICH SINE COLLAPSE WITH NGP INTERPOLATION # # # # # # 
    #box.rho[-1,:,:] = 1 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
        
    return box.rho
      
        

# --------------------------------------------- GREEN'S FUNCTION ---------------------------------------- #

def Greens_function(kx, ky, kz, a):
    """
    Returns Green's functions for discretised Poisson's equation
    """
    global Omega_0
    
    G = sum([(np.sin(i/2))**2 for i in [kx,ky,kz]])    
    G = -(3 * Omega_0)/(8 * a * G)
    
    return G

# ----------------------------------------- SOLVE POISSON EQUATION -------------------------------------- #

def get_phi(box, a):
    """
    Determine overdensity delta from rho. Then Fast Fourier transform (FFT) and multiply by corresponding Greens
    function. Inverse FFT to then obtain gravitational potential phi
    """
    
    # calculate overdesity delta
    delta = box.rho - np.ones((box.N_g, box.N_g, box.N_g))
    
    # FFT delta 
    FFT_delta = fft.rfftn(delta)
    
    # Test: inverse FFT of FFT of delta should be delta itself (i.e. IFFT[FFT(a)] = a) 
    # delta_test = fft.irfftn(fft.rfftn(delta), delta.shape)

    
    kx = 2 * np.pi * np.fft.fftfreq(delta.shape[0]) 
    ky = 2 * np.pi * np.fft.fftfreq(delta.shape[1]) 
    kz = 2 * np.pi * np.fft.rfftfreq(delta.shape[2])
    
    G = np.zeros(FFT_delta.shape)
    
    # cycle through each element of FFT_delta, multiply by Greens function, then inverse FFT to find 
    # gravitational potential field phi
    
    for l in range(len(kx)):
        for m in range(len(ky)):
            for n in range(len(kz)):

                kx_i = kx[l]
                ky_i = ky[m]
                kz_i = kz[n]
                
                # set phi[0,0,0] = 0
                if [kx_i,ky_i,kz_i] == [0,0,0]:
                    continue
                else:
                    G[l,m,n] = Greens_function(kx_i,ky_i,kz_i,a)
        
    
    FFT_phi = FFT_delta * G
    
    # Inverse FFT to obtain grav. potential phi        
    box.phi = np.fft.irfftn(FFT_phi)
    box.phi = box.phi

    return box.phi

# ------------------------------------- FORWARD INTEGRATE USING LEAPROG --------------------------------- #

def NGP_forward_integrate(box, a, da):
    
    f_a = (((1/a)*(box.Omega_M + box.Omega_K*a + box.Omega_Lambda*a*a*a))**-0.5)
    
    a_dash = a + (da / 2)
    f_a_dash = (((1/a_dash)*(box.Omega_M + box.Omega_K*a + box.Omega_Lambda*a*a*a))**-0.5)
    
    for i in range(len(box.x)):
        
        x_cell = int(box.x[i] // 1)
        y_cell = int(box.y[i] // 1)
        z_cell = int(box.z[i] // 1)
        
        if x_cell >= box.N_g:
            x_cell = x_cell % box.N_g    
        if x_cell == box.N_g - 1:
            x_cell_f = 0
        else:
            x_cell_f = x_cell + 1   

        if y_cell >= box.N_g:
            y_cell = y_cell % box.N_g 
        if y_cell == box.N_g - 1:
            y_cell_f = 0
        else:
            y_cell_f = y_cell + 1
            
        if z_cell >= box.N_g:
            z_cell = z_cell % box.N_g
        if z_cell == box.N_g - 1:
            z_cell_f = 0
        else:
            z_cell_f = z_cell + 1
        
        
        box.px[i] = box.px[i] - f_a * ((box.phi[x_cell_f, y_cell, z_cell] - box.phi[x_cell-1, y_cell, z_cell])/2) * da
        box.py[i] = box.py[i] - f_a * ((box.phi[x_cell, y_cell_f, z_cell] - box.phi[x_cell, y_cell-1, z_cell])/2) * da
        box.pz[i] = box.pz[i] - f_a * ((box.phi[x_cell, y_cell, z_cell_f] - box.phi[x_cell, y_cell, z_cell-1])/2) * da
        
        box.x[i] = box.x[i] + ((f_a_dash * box.px[i] * da) / (a_dash * a_dash))
        box.y[i] = box.y[i] + ((f_a_dash * box.py[i] * da) / (a_dash * a_dash))
        box.z[i] = box.z[i] + ((f_a_dash * box.pz[i] * da) / (a_dash * a_dash))
        
        if box.x[i] >= box.N_g:
            box.x[i] = box.x[i] % box.N_g
        if box.x[i] < 0:
            box.x[i] = box.N_g - (np.abs(box.x[i]) % box.N_g)
            
        if box.y[i] >= box.N_g:
            box.y[i] = box.y[i] % box.N_g
        if box.y[i] < 0:
            box.y[i] = box.N_g - (np.abs(box.y[i]) % box.N_g)
            
        if box.z[i] >= box.N_g:
            box.z[i] = box.z[i] % box.N_g
        if box.z[i] < 0:
            box.z[i] = box.N_g - (np.abs(box.z[i]) % box.N_g)
            
        
    return  

#I/O routines

def save_universe_details(box,address):
    
    file = open(address, "w") 
    file.write("------------ PM Simulation -------------\n") 
    file.write("Omega_0 = 1.0\n")
    file.write("Omega_M = {0}\n".format(box.Omega_M)) 
    file.write("Omega_K = {0}\n".format(box.Omega_K))
    file.write("Omega_Lambda = {0}\n".format(box.Omega_Lambda))
    file.write("-------------------------------------------------\n")
    file.write("Universe side length = {0} Mpc\n".format(box.L_box)) 
    file.write("Universe grid number = {0}".format(box.N_g)) 
    file.close()

def IC_output(Sx,Sy,Sz,Px,Py,Pz,address):
    
    f = h5py.File(address, 'w')
    f.create_dataset('Sx', data=Sx)
    f.create_dataset('Sy', data=Sy)
    f.create_dataset('Sz', data=Sz)
    f.create_dataset('Px', data=Px)
    f.create_dataset('Py', data=Py)
    f.create_dataset('Pz', data=Pz)
    f.close

def VTK_output_x_and_p(x,y,z,px,py,pz,address):
    
    coord = np.array([[x[i],y[i],z[i]] for i in range(len(x))])
    coord = list(coord.flatten())
    vels = np.array([[px[i],py[i],pz[i]] for i in range(len(px))])
    vels = list(vels.flatten())
    var = (("coord", 3, 1, coord), ("vels", 3, 1, vels))
    visit_writer.WritePointMesh(address, 0, coord, var)
    print ' VTK OUTPUT COMPLETE'
    
def VTK_output_rho(rho,N_g,address):
    
    nodal = list(rho.flatten())
    zonal = np.ndarray((N_g-1,N_g-1,N_g-1))
    for x in range(N_g-1):
        for y in range(N_g-1):
            for z in range(N_g-1):
                zonal[x,y,z] = (rho[x,y,z] + rho[x+1,y+1,z+1]) / 2
    zonal = list(zonal.flatten())
    dims = (N_g, N_g, N_g)
    var = (("zonal", 1, 0, zonal), ("nodal", 1, 1, nodal))
    visit_writer.WriteRegularMesh(address, 0, dims, var)
    
def HDF5_output_rho(rho,address):
    f = h5py.File(address, 'w')
    f.create_dataset('rho', data=rho)
    f.close
    
def HDF5_output_x_and_p(x,y,z,px,py,pz,address):
    
    f = h5py.File(address, 'w')
    f.create_dataset('x', data=x)
    f.create_dataset('y', data=y)
    f.create_dataset('z', data=z)
    f.create_dataset('px', data=px)
    f.create_dataset('py', data=py)
    f.create_dataset('pz', data=pz)
    f.close
    print 'output complete to hdf5'

    
    
    
    