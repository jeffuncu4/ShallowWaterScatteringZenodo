#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 19:04:08 2021

@author: jeff
"""


import numpy as np
import os
import subprocess
import pickle
from AnalysisClass import Analysis
from PlottingClass import Plot


class Simulation:
    def __init__(self, Ro, Bu, Lr, Ur, cyclonic = False):
        ###PARAMETERS are intialized as floats 
        
        #non dimensional numbers
        self.Ro = Ro
        self.Bu = Bu
        self.Lr = Lr
        self.Ur = Ur
        self.cyclonic = cyclonic
        
        ### SIMULATION LOGISITICS
        if self.cyclonic:
            self.vortex_name = 'Ro{}Bu{}_vortex_cyclonic.h5'.format(self.Ro, self.Bu)
            self.exp_name = 'Ro{}Bu{}Lr{}Ur{}_cyclonic'.format(self.Ro, self.Bu, self.Lr, self.Ur)
            
        else: 
            self.vortex_name = 'Ro{}Bu{}_vortex.h5'.format(self.Ro, self.Bu)
            self.exp_name = 'Ro{}Bu{}Lr{}Ur{}'.format(self.Ro, self.Bu, self.Lr, self.Ur)
        
        
        self.home_dir = '/home/jeff/Dedalus_projects/ITgeosims/ShallowWaterScatteringZenodo/'   #CHANGE THIS
        #where codes and bash script will be stored and where the git repo is
        
        self.bash_script_dir = self.home_dir + 'bash_scripts/'
        self.code_dir = self.home_dir + 'codes/'
        self.data_dir = '/media/jeff/Data/SWportableData/'  #CHANGE THIS
        #Where the data for each experiment is stored
        
        self.exps_dir = self.data_dir + 'experiments/'
        #where all the experiments are held
        self.vortex_dir = self.data_dir + 'vortices/'
        
        self.exp_dir = self.exps_dir + self.exp_name + '/'
        #where the specific directory is held
                        
        print (self.exp_name + ' initialized')                                          


        self.analysis = None
#        self.plots = Plot(Ro, Bu, Lr, Ur)
        self.plots = None
        

        
        self.Re = 1e6
        self.sponge_decay = 0.05
        self.forcing_decay = 1. #0.2
        
        #CONSTANTS
        self.g = 9.81 # m/s**2 
        
        #Dimensional numbers I choose
        self.f = 1e-4 # 1/s
        self.L = 25000 #m
        
        #Resulting dimensional quantities DEFUNE HG FIRST NOT SURE IF THIS I SMART
        self.wavelength = self.L/self.Lr
        self.hg = self.Ro*(self.L*self.f)**2/self.g
        self.ug  = self.hg/(self.L)*self.g/self.f # Not explicitly used in cyclogestorphic vortices
        self.H = self.Bu*(self.f*self.L)**2/self.g
        self.uw = self.ug/self.Ur
 
        self.l = np.pi*2/self.wavelength
        self.omega = np.sqrt(self.g*self.l**2*self.H + self.f**2)
        self.T_wave = np.pi*2/self.omega   
       
        self.hw = self.l*self.H*self.uw/self.f
        self.vw = self.uw*self.omega/self.f
        self.tau_s = self.sponge_decay*2*np.pi/self.omega
        self.tau_w = self.forcing_decay*2*np.pi/self.omega
        self.c = np.sqrt(self.g*self.H)
        self.Ld = self.c/self.f
      
        
        
        #Simulations  non dimensional numbers
        self.L_aspect = 1
        self.forcing_window_length_wavelength = 1.
        self.sponge_window_length_wavelength = 1.
        self.geo_window_length_ratio = 0.2 # doesnt do anyting rn
        
#        if self.Lr >=1:
#            self.nx = 640
#        else:
#            self.nx = 1280
        self.nx= 512   
        
#        if self.Lr >=1:
#            self.nx = 512
#        else:
#            self.nx = 1024
        # I choose to define dx, dt for the simulation Nondim# 
        #Keep in mind I neeed to satisfy some sort of cfl condition with dt <  C*dx/U
        self.CFL_constant = 0.002
        
        self.dx = 500 # m
        self.dt_cfl = self.CFL_constant*self.dx/self.ug # some sort of round function to nearest integer divsion of wave period
        # because wave_period == dt*N whre is N is integer must be true for flux calc
         
            
        #
        
        def round_dt(dt_orig, wave_period):
            ''' 
            This function allows for wave_period./dt=integer. T
            his is useful for Fourier Transforms
            '''
            N = 8 # factor by which all periods will be divisible by
            save_iter = 0
            dt_diff = 1.
            new_dt = None
            while dt_diff >= 0 :
                wave_fraction = wave_period/N
                #print (N, wave_fraction)
                dt_diff = wave_fraction - dt_orig
                new_dt = wave_fraction
                N+= 8
                save_iter+= 1
                #print (save_iter)
            return new_dt, save_iter
        
        
        self.dt, self.save_iter  = round_dt(self.dt_cfl, self.T_wave)
    
        
        # The resulting parameters are
        self.mu = self.ug*self.dx**3/self.Re  #m^4/s
        # not true if you use CFL
    
        
        self.Lx = self.nx*self.dx
        self.Ly = self.Lx*self.L_aspect
        self.ny = self.nx*self.L_aspect
        self.forcing_window_length = self.forcing_window_length_wavelength*(np.pi*2/self.l)
        self.sponge_window_length = self.sponge_window_length_wavelength*(np.pi*2/self.l)
        self.geo_window_length = self.geo_window_length_ratio*(np.pi*2/self.l) # dont think I m using this
        
        
        self.inertial_period = np.pi*2/self.f
        self.c_rot = self.omega/self.l
        self.num_periods = self.Ly*4./3./self.c_rot/self.inertial_period 
        self.Ts = self.Ly/self.c_rot*4./3
        
        #Other Numbers which are relevant 
             
        

        
        #self.save_iter = 10 # function of wave period main sim save iteration
        self.num_iter = int(self.Ts/self.dt)
        

        self.num_iter_vortex = int(self.inertial_period*3/self.dt)
        self.save_iter_vortex= self.num_iter_vortex//4# I really only need the last frame
        
        self.max_writes = 100000        
        
        self.nx_v = 512
        self.Lx_v = self.nx_v*self.dx
        self.Ly_v = self.Lx_v*self.L_aspect
        self.ny_v = self.nx_v*self.L_aspect
        
        

    def create_sim(self, run_if_created = False):       
        sim_created = self.exp_name in os.listdir(self.exps_dir)
        if sim_created and not run_if_created:
            print ("Simulation already created")
            return None
        
        subprocess.call(self.bash_script_dir +"create_sim.sh " + self.exp_name + 
                        ' ' + self.exps_dir,  shell=True) # add arguement for data dir
        print ("creating simulation")
            
        
        parameterValues = [self.Ur,self.Lr,self.Ro , self.Bu, self.Re, 
                           self.sponge_decay, self.forcing_decay,self.ug ,self.f,self.g,self.uw, self.L, self.l,
                self.H ,self.mu, self.omega, self.hg ,self.hw ,self.vw ,self.tau_s,self.tau_w ,self.num_periods,
                self.L_aspect,self.forcing_window_length_wavelength,self.sponge_window_length_wavelength ,
               self. geo_window_length_ratio ,self.dt,self.dx, self.Ts ,self.num_iter, self.Lx ,self.Ly ,self.nx,
                self.forcing_window_length,self.sponge_window_length,self.geo_window_length ,
                self.save_iter ,self.ny ,self.T_wave,self.wavelength,self.c ,self.Ld ,self.max_writes, self.num_iter_vortex, 
                self.save_iter_vortex, self.nx_v, self.Lx_v, self.Ly_v, self.ny_v, self.exp_dir, self.vortex_name,
                self.cyclonic]


        parameterNames = ['Ur','Lr','Ro' ,'Bu' ,'Re','sponge_decay','forcing_decay','ug' ,'f','g','uw', 'L', 'l',
                'H' ,'mu', 'omega','hg' ,'hw' ,'vw' ,'tau_s','tau_w' ,'num_periods',
                'L_aspect','forcing_window_length_wavelength','sponge_window_length_wavelength' ,
                'geo_window_length_ratio' ,'dt','dx', 'Ts' ,'num_iter', 'Lx' ,'Ly' ,'nx',
                'forcing_window_length','sponge_window_length','geo_window_length' ,
                'save_iter' ,'ny' ,'T_wave','wavelength','c' ,'Ld' ,'max_writes',
                'num_iter_vortex', 'save_iter_vortex',  'nx_v', 'Lx_v', 'Ly_v', 'ny_v', 'exp_dir', 'vortex_name', 
                'cyclonic']
        
        my_dict = dict(zip(parameterNames, parameterValues))
        with open(self.exp_dir + 'IC/parameters.pkl', 
                  'wb') as pfile:
            pickle.dump(my_dict, pfile)
            
        
        subprocess.call('python forcing_windows.py ' + self.exp_dir,  shell=True) # could add extra argument here too

   
    def run_vortex_sim(self, run_if_created = False):
        ''' this will run the cylogeostophic adjustment code as and run the simulation
        until it settles'''
        
        
        vortex_created  = self.vortex_name in os.listdir(self.vortex_dir)
        
        #vortex_created = 'settled_vortex.h5' in os.listdir(self.exp_dir + 'IC/')
        if vortex_created and not run_if_created:
            print ("vortex already created already created")
            subprocess.call('cp ' + self.vortex_dir + self.vortex_name + ' '
                            + self.exp_dir + 'IC/' + self.vortex_name, shell=True)
            subprocess.call('python array_inception.py ' + self.exp_dir , shell=True)
            print (' copied vortex into experiment folder')
            return None
        
        subprocess.call(self.bash_script_dir + "run_vortex_sim.sh {}".format(self.exp_dir), shell=True)
    
        
    def run_main_sim(self, run_if_created = False): #this creates and runs the simulation
        # this should check if the siulation already exist in vortex_fields
        data_created = 'data_s1.h5' in os.listdir(self.exp_dir + 'data/')
        if data_created and not run_if_created:
            print ("main sim has been run")
            self.analysis = Analysis(self.Ro, self.Bu, self.Lr, self.Ur, cyclonic = self.cyclonic)
            return None
        subprocess.call(self.bash_script_dir + "run_main_sim.sh {}".format(self.exp_dir), shell=True)
        self.analysis = Analysis(self.Ro, self.Bu, self.Lr, self.Ur, cyclonic = self.cyclonic)
    

    def make_uvh_videos(self, run_if_created=False):
        vids_created = 'h.mp4' in os.listdir(self.exp_dir + 'plots/')
        if vids_created and not run_if_created:
            print ("uvh videos have been created")
            return None
        subprocess.call(self.bash_script_dir + "make_uvh_videos.sh {}".format(self.exp_dir), shell=True)
    
    def run_analysis(self):
        pass

    def run_sim(self, run_if_created=False):
        self.create_sim(run_if_created)
        self.run_vortex_sim(run_if_created)
        self.run_main_sim(run_if_created)
        #self.make_uvh_videos(run_if_created)


#
#
#    def run_and_analyze(self):
#        self.run_all()
#        self.make_videos()
#        self.run_analysis()


#    def run_analytics(self):
#        if self.sim_executed == True:
#            subprocess.call("run_analytics.sh {}".format(self.expName))
#            
#    def plot_all(self):
#        if self.sim_executed == True:
#            subprocess.call("plot_all.sh {}".format(self.expName))
#        
