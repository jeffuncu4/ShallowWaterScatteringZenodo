#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 19:48:23 2021

@author: jeff
"""
import pickle
import numpy as np
import h5py
import numpy.fft as fft
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import cm
#import scipy.fftpack as sfft
#
class Analysis:
    def __init__(self, Ro, Bu, Lr, Ur, cyclonic = False):
    ###PARAMETERS
    
        #non dimensional numbers
        self.Ro = Ro
        self.Bu  = Bu
        self.Lr = Lr
        self.Ur = Ur
        self.cyclonic = cyclonic
        

        if self.cyclonic:
            self.vortex_name = 'Ro{}Bu{}_vortex_cyclonic.h5'.format(self.Ro, self.Bu)
            self.exp_name = 'Ro{}Bu{}Lr{}Ur{}_cyclonic'.format(self.Ro, self.Bu, self.Lr, self.Ur)
            
        else: 
            self.vortex_name = 'Ro{}Bu{}_vortex.h5'.format(self.Ro, self.Bu)
            self.exp_name = 'Ro{}Bu{}Lr{}Ur{}'.format(self.Ro, self.Bu, self.Lr, self.Ur)
        
        self.home_dir = '/home/jeff/Dedalus_projects/ITgeosims/ShallowWaterScatteringZenodo/' #CHANGE THIS
        self.bash_script_dir = self.home_dir + 'bash_scripts/'
        self.code_dir = self.home_dir + 'codes/'
        self.data_dir = '/media/jeff/Data/SWportableData/'  #CHANGE THIS
        self.exps_dir = self.data_dir + 'experiments/'
       

        self.exp_dir = self.exps_dir + self.exp_name + '/'
        self.data_file_name = self.exp_dir + 'data/data_s1.h5'
        self.vortex_dir = self.data_dir + 'vortices/'


        self.par_dict = pickle.load(open(self.exp_dir + 'IC/parameters.pkl', 'rb'))
        
        self.f = self.par_dict['f']
        self.L = self.par_dict['L']
        
        
        self.g = self.par_dict['g']
        self.Ly = self.par_dict['Ly']

        self.omega = self.par_dict['omega']
        self.l = self.par_dict['l']
        self.wavelength = self.par_dict['wavelength']
        self.H = self.par_dict['H']
        
        
        self.c_rot = self.omega/self.l
        
        data_file = h5py.File(self.data_file_name, mode = 'r')
        self.u = np.array(data_file['tasks']['u'])
        self.v = np.array(data_file['tasks']['v'])
        self.h = np.array(data_file['tasks']['h'])
        self.time_axis = np.array(data_file['scales']['sim_time'])
        self.x_axis = np.array(data_file['scales']['x']['1.0'])
        self.y_axis = np.array(data_file['scales']['y']['1.0'])
        data_file.close() 
        self.dt = self.time_axis[1]-self.time_axis[0] # note this is not the dt the simulation ran at
        self.dx = self.x_axis[1] - self.x_axis[0]
        self.nt = len(self.time_axis)
        self.nx = len(self.x_axis)
        self.ny = len(self.y_axis)
        self.shape = np.shape(self.u)
        
        
        self.k_array = fft.fftshift(fft.fftfreq(self.nx, d=self.dx))*np.pi*2
        self.l_array = fft.fftshift(fft.fftfreq(self.ny, d=self.dx))*np.pi*2
        self.omega_array = fft.fftshift(fft.fftfreq(self.nt, d=self.dt)*2 *np.pi)
        print (np.shape(self.u))
        
        
        self.u0 = self.u[0]
        self.v0 = self.v[0]
        self.h0 = self.h[0]
        
        
        self.h_wave = self.h-self.h0
        self.u_wave = self.u-self.u0
        self.v_wave = self.v-self.v0
    
    def vorticity(self):
        vx = np.gradient(self.v[0], axis = 0)/self.dx
        uy = np.gradient(self.u[0], axis = 1)/self.dx
        return vx - uy
    

        
        
        
        
        
        
        
        