#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 09:58:02 2021

@author: jeff
"""

from simulationClass import Simulation 
import numpy as np
import matplotlib.pyplot as plt
from AnalysisClass import Analysis



Ro = 0.04 # must make these floats for proper naming conventions
Bu = 1.0
Lr = 3.0
Ur = 1000.0
cyclonic=True

def adjusted_params(Ro, Bu, Lr, cyclonic):
    ana = Analysis(Ro, Bu, Lr, 1000.0, cyclonic=cyclonic)
    print (ana.exp_name)
    exp_name  = ana.exp_name
    L = ana.L
    x = ana.x_axis
    f = ana.f
    u_cross = ana.u[0, 256, :]
    u_max = np.max(u_cross)
    umi = np.argmax(u_cross)
    R = np.abs(x[umi])
    print (f'scaled ratio {R/(L/np.pi)}')
    L_adj = np.pi*R
    lam_orig = L/Lr
    K_adj = L_adj/lam_orig
    # Bu_adj = Bu/L**2*L_adj**2
    Bu_adj = Bu*L**2/L_adj**2
    # Bu_adj = Bu*1
    Ro_bulk_adj = u_max/L_adj/f
    print (f'k_adj: {K_adj}, Bu_adj: {Bu_adj}, Ro: {Ro}')
    local_Ro = np.max(np.abs(ana.vorticity()))/ana.f
    print( exp_name, Ro_bulk_adj, local_Ro, Bu_adj, K_adj, L_adj/L)

adjusted_params(Ro, Bu, Lr, cyclonic)




        