"""
Created on Wed Jan 10 15:42:34 2018

@author: James Rogers
"""

from __future__ import division
import numpy as np
from numpy import pi,sin,cos
import scipy as sp

ns = 0.96                  # Spectral index
H0_normal_units = (70)     # Hubble constant in c =/= 1 units

def window(x):
    
    W = (sin(x) - x*cos(x)) / (x**3)
    
    return 3*W

def sigma_integrand(k):
    
    W = window((800/H0_normal_units) * k)
    I = (W**2) * (k**ns) * (k**2)
    return I


def sigma_0(kmax):
    sigma_0 = sp.integrate.quad(sigma_integrand, 1e-7, kmax)
    return ((4*pi))*sigma_0[0]

###############################
#N_g = 128
#r_0 = 1000/N_g
#
#kmax = 2 * pi
#sigma0 = sigma_0(kmax)
#print 'sigma0 =', sigma0

#print 'As =', 0.001 * (0.9*0.9) / sigma0