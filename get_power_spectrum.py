"""
Created on Fri Dec 01 13:42:47 2017

@author: James Rogers
"""

import h5py
import numpy as np

N_g = 32
#simulation = "OUTPUT\04.12.2017_20.10.39"

f = h5py.File("collated_data", "r")
rho = f["rho_final"][:]

delta = rho - np.ones((N_g, N_g, N_g))
FFT_delta = np.fft.rfftn(delta)
FFT_delta_squared = FFT_delta*FFT_delta

kx = 2 * np.pi * np.fft.fftfreq(delta.shape[0]) 
ky = 2 * np.pi * np.fft.fftfreq(delta.shape[1]) 
kz = 2 * np.pi * np.fft.rfftfreq(delta.shape[2]) 

k_data = []

for x in range(len(kx)):
    for y in range(len(ky)):
        for z in range(len(kz)):
            
            k = np.sqrt((kx[x]*kx[x] + ky[y]*ky[y] + kz[z]*kz[z]))
            k_values = [k_data[i][0] for i in range(len(k_data))]
            
            if k not in k_values:
                k_data.append([k, [FFT_delta_squared[x,y,z]]])
            else:
                index = k_values.index(k)
                k_data[index][1].append(FFT_delta_squared[x,y,z])
                
print len(k_data)
print k_data[3]