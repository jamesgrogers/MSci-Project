# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 17:29:20 2017

@author: James Rogers
"""

from __future__ import division
import numpy as np


def CIC_density_assignment(box):
    """
    Cloud in Cell density assignment assumes particles box shaped and particle's mass is 
    assigned across neighbouring cells.
    """
    N_g = box.N_g
    box.rho = np.zeros((box.N_g, box.N_g, box.N_g))
    
    for i in range(len(box.x)):
        x_cell = int(box.x[i] // 1) % N_g
        y_cell = int(box.y[i] // 1) % N_g
        z_cell = int(box.z[i] // 1) % N_g
        
        #print (x_cell,y_cell,z_cell)
        
        x_cell_f = (x_cell + 1) % N_g
        y_cell_f = (y_cell + 1) % N_g
        z_cell_f = (z_cell + 1) % N_g
        
        #print (x_cell_f,y_cell_f,z_cell_f)
        
        dx = (box.x[i] - x_cell) % N_g
        dy = (box.y[i] - y_cell) % N_g
        dz = (box.z[i] - z_cell) % N_g
        
        #print (dx,dy,dz)
        
        tx = 1 - dx
        ty = 1 - dy
        tz = 1 - dz
        
        #print (tx,ty,tz)
        
        box.rho[x_cell, y_cell, z_cell] += box.m[i]*tx*ty*tz
        
        box.rho[x_cell_f, y_cell, z_cell] += box.m[i]*dx*ty*tz
        box.rho[x_cell, y_cell_f, z_cell] += box.m[i]*tx*dy*tz
        box.rho[x_cell, y_cell, z_cell_f] += box.m[i]*tx*ty*dz
        
        box.rho[x_cell_f, y_cell_f, z_cell] += box.m[i]*dx*dy*tz
        box.rho[x_cell_f, y_cell, z_cell_f] += box.m[i]*dx*ty*dz
        box.rho[x_cell, y_cell_f, z_cell_f] += box.m[i]*tx*dy*dz 
        
        box.rho[x_cell_f, y_cell_f, z_cell_f] += box.m[i]*dx*dy*dz
        
############################## FORWARD INTEGRATE ###############################      

def phi_gradient_x(box,i,j,k):
    
    i_f = (i+1) % box.N_g
    gx = -(box.phi[i_f,j,k] - box.phi[i-1,j,k])/2
    return gx

def phi_gradient_y(box,i,j,k):
    
    j_f = (j+1) % box.N_g
    gy = -(box.phi[i,j_f,k] - box.phi[i,j-1,k])/2
    return gy

def phi_gradient_z(box,i,j,k):
    
    k_f = (k+1) % box.N_g
    gz = -(box.phi[i,j,k_f] - box.phi[i,j,k-1])/2
    return gz

def CIC_gradient(box,
                 x_cell, y_cell, z_cell,
                 x_cell_f, y_cell_f, z_cell_f,
                 tx, ty, tz, dx, dy, dz):
    
    gp_x = phi_gradient_x(box,x_cell,y_cell,z_cell) * tx*ty*tz\
         + phi_gradient_x(box,x_cell_f,y_cell,z_cell) * dx*ty*tz\
         + phi_gradient_x(box,x_cell,y_cell_f,z_cell) * tx*dy*tz\
         + phi_gradient_x(box,x_cell,y_cell,z_cell_f) * tx*ty*dz\
         + phi_gradient_x(box,x_cell_f,y_cell_f,z_cell) * dx*dy*tz\
         + phi_gradient_x(box,x_cell_f,y_cell,z_cell_f) * dx*ty*dz\
         + phi_gradient_x(box,x_cell,y_cell_f,z_cell_f) * tx*dy*dz\
         + phi_gradient_x(box,x_cell_f,y_cell_f,z_cell_f) * dx*dy*dz
         
    gp_y = phi_gradient_y(box,x_cell,y_cell,z_cell) * tx*ty*tz\
         + phi_gradient_y(box,x_cell_f,y_cell,z_cell) * dx*ty*tz\
         + phi_gradient_y(box,x_cell,y_cell_f,z_cell) * tx*dy*tz\
         + phi_gradient_y(box,x_cell,y_cell,z_cell_f) * tx*ty*dz\
         + phi_gradient_y(box,x_cell_f,y_cell_f,z_cell) * dx*dy*tz\
         + phi_gradient_y(box,x_cell_f,y_cell,z_cell_f) * dx*ty*dz\
         + phi_gradient_y(box,x_cell,y_cell_f,z_cell_f) * tx*dy*dz\
         + phi_gradient_y(box,x_cell_f,y_cell_f,z_cell_f) * dx*dy*dz
         
    gp_z = phi_gradient_z(box,x_cell,y_cell,z_cell) * tx*ty*tz\
         + phi_gradient_z(box,x_cell_f,y_cell,z_cell) * dx*ty*tz\
         + phi_gradient_z(box,x_cell,y_cell_f,z_cell) * tx*dy*tz\
         + phi_gradient_z(box,x_cell,y_cell,z_cell_f) * tx*ty*dz\
         + phi_gradient_z(box,x_cell_f,y_cell_f,z_cell) * dx*dy*tz\
         + phi_gradient_z(box,x_cell_f,y_cell,z_cell_f) * dx*ty*dz\
         + phi_gradient_z(box,x_cell,y_cell_f,z_cell_f) * tx*dy*dz\
         + phi_gradient_z(box,x_cell_f,y_cell_f,z_cell_f) * dx*dy*dz
    
    return gp_x,gp_y,gp_z

    
def CIC_forward_integrate(box, a, da):

    N_g = box.N_g
    
    f_a = (((1/a)*(box.Omega_M + box.Omega_K*a + box.Omega_Lambda*a*a*a))**-0.5)
    
    a_dash = a + (da / 2)
    f_a_dash = (((1/a_dash)*(box.Omega_M + box.Omega_K*a + box.Omega_Lambda*a*a*a))**-0.5)
    
    for i in range(len(box.x)):
        x_cell = int(box.x[i] // 1) % N_g
        y_cell = int(box.y[i] // 1) % N_g
        z_cell = int(box.z[i] // 1) % N_g
        
        x_cell_f = (x_cell + 1) % N_g
        y_cell_f = (y_cell + 1) % N_g
        z_cell_f = (z_cell + 1) % N_g
        
        dx = (box.x[i] - x_cell) % N_g
        dy = (box.y[i] - y_cell) % N_g
        dz = (box.z[i] - z_cell) % N_g
        
        tx = 1 - dx
        ty = 1 - dy
        tz = 1 - dz
        
        gx,gy,gz = CIC_gradient(box,x_cell,y_cell,z_cell,x_cell_f,y_cell_f,z_cell_f,tx,ty,tz,dx,dy,dz)
        
        box.px[i] = box.px[i] + f_a * gx * da
        box.py[i] = box.py[i] + f_a * gy * da
        box.pz[i] = box.pz[i] + f_a * gz * da
        
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