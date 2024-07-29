#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 11:14:39 2022

@author: jeff
"""

from simulationClass import Simulation 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import itertools as it
from matplotlib import cm
from matplotlib.colors import TwoSlopeNorm
from matplotlib.colors import Normalize
from AnalysisClass import Analysis



plt.rcParams.update({  # LaTeX fonts
    "text.usetex": True,
    "font.family": "Helvetica",
    "font.size": 12,
    'text.latex.preamble': r"\usepackage{amsmath}"
})

ro = 0.03
bu = 1.0
lr = 4.0
ur = 1000.0
cyclonic =False
exp = Simulation(ro, bu, lr, ur)
exp.run_sim()

Lx = exp.Lx
Ly = exp.Ly
Ld = exp.Ld
f = exp.f
L = exp.L




ana = Analysis(ro, bu, lr, 1000.0, cyclonic=cyclonic)
L = ana.L
x = ana.x_axis
f = ana.f
u_cross = ana.u[0, 256, :]
u_max = np.max(u_cross)
umi = np.argmax(u_cross)
R = x[umi]
R_scaled = R/L

print (ana.c_rot)

cg = ana.g*ana.H*ana.l/ana.omega


print (f'cg = {cg}, cp = {ana.c_rot}')

ana = exp.analysis
vort  = -1*ana.vorticity()/f #analyisclass has vort multiplied by negative one


vmin = np.min(vort)
vmax = np.max(vort)
norm = TwoSlopeNorm(vmin=vmin*2.5, vcenter=0, vmax=vmax)

u = ana.v0 #in analysis u is switched with v and f goes to soutpole 
v = ana.u0

extent = np.array([-Ly/2, Ly/2, -Lx/2, Lx/2])/L
print (1.25/extent*0.5)


fig, ax = plt.subplots(1, 1)

plt.imshow(vort, extent = extent, cmap = 'PRGn', norm = norm, origin='lower')
plt.ylabel(r'$y/L$', fontsize = 'x-large')
plt.xlabel(r'$x/L$', fontsize = 'x-large')
clb = plt.colorbar()
clb.ax.set_title(r'$\zeta/f$', fontsize = 'x-large')
# plt.show()

# us = u[ana.nx//2, :]

# print (np.where(us == np.max(us)))
# plt.imshow(u, extent = extent)
# plt.plot([1/np.pi, 1/np.pi], [-5, 5], 'k--')
# plt.show()



slicevel = slice(0, len(u), 4)
x = np.arange(-64, 64)/128*Lx/L
ucrop = u[slicevel, slicevel]
vcrop = v[slicevel, slicevel]

colors = (ucrop**2 + vcrop**2)**0.5

norm = Normalize()
norm.autoscale(colors)
colormap = cm.copper


plt.quiver(x, x,ucrop, vcrop, scale=1, scale_units='inches')
# plt.plot([-0.5, -0.5], [-5, 0], 'k--')
# plt.plot([0.5, 0.5], [-5, 0], 'k--')
plt.plot([0,0], [0, 5], 'k--')
plt.plot([R_scaled, R_scaled], [0, 5], 'k--')



plt.xlim(-1.25, 1.25)
plt.ylim(-1.25, 1.25)


# ax.annotate("", xy=(-0.5, -0.8), xytext=(0.5, -0.8),
#             arrowprops=dict(arrowstyle="<->"))
# ax.text(0, -1, "$L$", horizontalalignment='center',
#         verticalalignment='bottom', bbox=dict(boxstyle="Round", fc="w", alpha=0.5, lw=0.) , fontsize = 'x-large')


ax.annotate("", xy=(0, 1), xytext=(R_scaled, 1),
            arrowprops=dict(arrowstyle="<->"))
ax.text(0.11, 0.75, "$R$", horizontalalignment='center',
        verticalalignment='bottom', bbox=dict(boxstyle="Round", fc="w", alpha=0.5, lw=0.) , fontsize = 'x-large')


plt.show()
