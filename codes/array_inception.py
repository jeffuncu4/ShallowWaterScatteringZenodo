#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 09:10:00 2021

@author: jeff
"""


import numpy as np
import h5py
import pickle
import sys
import matplotlib.pyplot as plt



def implant_vortex(vortex_data, nx_sim):
    #nx_sim must be even , len )vortex data must be even, and both must be square
    nx_vortex = len(vortex_data)
    expanded_vortex = np.zeros([nx_sim, nx_sim])
    expanded_vortex[nx_sim//2-nx_vortex//2: nx_sim//2 + nx_vortex//2, 
                    nx_sim//2-nx_vortex//2: nx_sim//2 + nx_vortex//2]  = vortex_data
    
    return expanded_vortex


exp_dir = sys.argv[1] 

par_dict = pickle.load(open(exp_dir + 'IC/parameters.pkl', 'rb'))



Ro = par_dict['Ro']
Bu = par_dict['Bu']
vortex_name = par_dict['vortex_name']

vortex_file = h5py.File(exp_dir + 'IC/' + vortex_name, 'r')



vortex_h = np.array(vortex_file.get('geoH'))
vortex_u = np.array(vortex_file.get('geoU'))
vortex_v = np.array(vortex_file.get('geoV'))
H = par_dict['H']
#print (H)

#
#plt.imshow(vortex_h-H)
#
#plt.title('h not expanded')
#plt.show()

nx = par_dict['nx']



expanded_vortex_h = implant_vortex(vortex_h-H, nx)  + H
expanded_vortex_u = implant_vortex(vortex_u, nx)
expanded_vortex_v = implant_vortex(vortex_v, nx)
#
#plt.imshow(expanded_vortex_h)
#plt.title('h_expanded')
#plt.savefig(exp_dir + 'IC/settled_vortex.png')
#plt.show()
#plt.imshow(expanded_vortex_u)
#plt.title('u_expanded')
#plt.show()
#



saveData=True
if saveData:
    hf = h5py.File(exp_dir + 'IC/' + 'expanded_' + vortex_name, 'w')
    hf.create_dataset('geoH', data=expanded_vortex_h)
    hf.create_dataset('geoU', data=expanded_vortex_u)
    hf.create_dataset('geoV', data=expanded_vortex_v)
    hf.close()



