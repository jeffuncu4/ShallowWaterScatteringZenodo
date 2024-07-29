#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 19:42:02 2021

@author: jeff
"""
import subprocess


class Plot:
    def __init__(self, Ro, Bu, Lr, Ur):
        self.Ro = Ro
        self.Bu = Bu
        self.Lr = Lr
        self.Ur = Ur
        self.bash_script_dir = '../bash_scripts/'
        ### SIMULATION LOGISITICS
        self.exp_name = 'Ro{}Bu{}Lr{}Ur{}'.format(self.Ro, self.Bu, self.Lr, self.Ur)
        self.exp_dir = '../experiments/' + self.exp_name + '/'
        self.analysis_dir = self.exp_dir + 'analysis/'
        self.plots_dir = self.exp_dir + 'plots/'
                                                                  
    
    def view(self, data_name, video = True):
        data_location = {'u':  'u.mp4',  'v': 'v.mp4',
                          'h': 'h.mp4'  }
        
        if video:
            subprocess.call('ffplay ' + self.plots_dir + data_location[data_name],  shell=True)
        else:
            subprocess.call('feh' + self.plots_dir + data_location[data_name],  shell=True)