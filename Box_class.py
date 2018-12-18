"""
MSci Project - PM code
@author: James Rogers

All dimensionless variables follow that of Anatoly Klypin's PM code 
"""
from __future__ import division
import numpy as np


class Box:
    def __init__(self, N_g, L_box, Omega_M, Omega_K, Omega_Lambda):
        """
        N_g should be number of cells in each axis e.g. N_g=64
        """
        # define box dimensions, side length and cell width r_0
        self.N_g = N_g
        self.L_box = L_box
        self.r_0 = L_box/N_g
        
        
        #define particle positions
        self.x = []
        self.y = []
        self.z = []
        
        #define particle velocities
        self.px = []
        self.py = []
        self.pz = []
        
        #define box density
        self.rho = np.zeros((self.N_g, self.N_g, self.N_g))
        
        #define gravitational potential phi
        self.phi = None
        
        #define cosmological parameters
        self.Omega_M = Omega_M
        self.Omega_K = Omega_K
        self.Omega_Lambda = Omega_Lambda
        
        # define particle masses
        self.m = [self.Omega_M]*(self.N_g)**3

        

    