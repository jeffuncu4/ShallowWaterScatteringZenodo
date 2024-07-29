#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 13:17:32 2022

@author: jeff
"""


from simulationClass import Simulation 
from AnalysisClass import Analysis
import numpy as np
import matplotlib.pyplot as plt


#Note these nunmbers will not be the final parameters after the cyclogeostrophic adjustment
Ro = 0.04 
Bu = 1.0
Lr = 3.0
Ur = 1000.0
cyclonic=True

#list of values used in the article
# ro_cyclonic_list = [0.01, 0.03, 0.04, 0.06] (cyclonic=True)
# ro_anticyclonic_list = [0.006, 0.01, 0.02, 0.03] (cyclonic=False)
# bu_list = [0.5, 0.9, 1.0, 1.1, 1.5]
# lr_list= [1.,1.5,  2., 3., 4.]

#Create and run simulation
exp = Simulation(Ro, Bu, Lr, Ur, cyclonic=cyclonic)
exp.run_sim()

#In case one wants to do indidivual parts of running this experiment
# exp.create_sim()
# exp.run_vortex_sim()
# exp.run_main_sim()

# Example for how to use Analysis class
ana = Analysis(Ro, Bu, Lr, Ur, cyclonic=cyclonic)

h_wave = ana.h- ana.h[0]
x = ana.x_axis
L = ana.L
extent = np.array([x[0], x[-1], x[0], x[-1]])/L

plt.imshow(h_wave[-1], extent = extent)
plt.title('Wave Height Field')
plt.xlabel('x/L')
plt.ylabel('y/L')
plt.show()