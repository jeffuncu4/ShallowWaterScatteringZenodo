#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 12:10:17 2020

@author: jeff
"""

import matplotlib.pyplot as plt
import numpy as np
import h5py
from dedalus import public as de
from scipy import signal
import sys
import pickle

#exp_dir = '../experiments/' + sys.argv[1] + '/'

exp_dir = sys.argv[1] 

par_dict = pickle.load(open(exp_dir + 'IC/parameters.pkl', 'rb'))

dx = par_dict['dx']
nx = par_dict['nx']
ny = par_dict['ny']
Lx = par_dict['Lx']
Ly = par_dict['Ly']
forcing_window_length = par_dict['forcing_window_length']
sponge_window_length = par_dict['sponge_window_length']



xbasis = de.Fourier('x', nx, interval=(-Lx/2,Lx/2), dealias=1)
ybasis = de.Fourier('y', ny, interval=(-Ly/2,Ly/2), dealias=1)
domain = de.Domain([xbasis, ybasis], grid_dtype=np.float64)
x = domain.grid(0)[0]
y = domain.grid(1)[0]

num_ppts_forcing = int (forcing_window_length/dx)
num_ppts_sponge = int (sponge_window_length/dx)



#def gaussianForcingWindow(shift, thickness):
#    return np.exp(-((y + shift*Ly/2)/(thickness*Ly/2.))**2)

def tukeyForcingWindow(tukey_length, shape):
    tw = signal.tukey(tukey_length, alpha = shape)
    zeropadding = np.zeros(ny-tukey_length)
    tw = np.append(tw, zeropadding)
    return tw

alpha = 0.7

waveForcingWindow = tukeyForcingWindow(num_ppts_forcing, alpha)
spongeWindow = tukeyForcingWindow(num_ppts_sponge, alpha)[::-1]



plot = False
if plot:
    plt.plot(y, waveForcingWindow, label="wave forcing")
    plt.plot(y, spongeWindow, label = 'sponge forcing')
    plt.legend()
    plt.show()

waveForcingWindow2d = np.zeros([nx, ny])
spongeWindow2d = np.zeros([nx, ny])

for i in range(nx):
    waveForcingWindow2d[i, :] = waveForcingWindow
    spongeWindow2d[i, :] = spongeWindow

wwf = domain.new_field(name = 'wave window')
wwf['g'] =waveForcingWindow2d



saveData=True
if saveData:
    hf = h5py.File(exp_dir + 'IC/forcingWindows.h5', 'w')
    hf.create_dataset('waveForcingWindow', data=waveForcingWindow2d)
    hf.create_dataset('spongeWindow', data=spongeWindow2d)
    hf.close()












