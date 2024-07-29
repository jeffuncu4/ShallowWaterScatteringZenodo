#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 13:11:31 2022

@author: jeff
"""


from simulationClass import Simulation 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import itertools as it
from matplotlib import cm
from matplotlib.colors import TwoSlopeNorm

Ro = 0.03 # must make these floats for proper naming conventions
Bu = 1.0
Lr = 4.0
Ur = 1000.0

exp = Simulation(Ro, Bu, Lr, Ur)
exp.run_sim()

Lx = exp.Lx
Ly = exp.Ly
Ld = exp.Ld
f = exp.f
H = exp.H
#hw = exp.hw

ana = exp.analysis

hw = 0.000218

h  = (ana.h[-1] - ana.h[0])/hw


vmin = np.min(h)
vmax = np.max(h)

print (vmin, vmax, H)
absmax = max(np.abs(vmin), vmax)
norm = TwoSlopeNorm(vmin= -absmax, vcenter=0, vmax= absmax)


extent = np.array([-Ly/2, Ly/2, -Lx/2, Lx/2])/Ld
plt.imshow(h, extent = extent, cmap = 'bwr', norm = norm)
plt.ylabel(r'$y/L_d$')
plt.xlabel(r'$x/L_d$')
clb = plt.colorbar()
clb.ax.set_title(r'$(h-H)/H_w$')

plt.show()
