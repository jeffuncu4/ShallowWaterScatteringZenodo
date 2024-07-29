from simulationClass import Simulation 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import itertools as it
from matplotlib import cm
import numpy.fft as fft
from AnalysisClass import Analysis
from scipy import interpolate
from matplotlib.colors import TwoSlopeNorm
import pandas as pd
from scipy.interpolate import RegularGridInterpolator as RGI

'''
Calculates the scattering ratio S for a given experiment with the final function "flux_ratio". 
This function is imported to save_experiment_statistics.py for calculating S for all experiments.
'''

def load_exp_vars(Ro, Bu, Lr, cyclonic,  Ur =1000.):
    '''output is a 3d varaiable for each field whose time axis goes
    from the nt- N*T to nt
    '''

    ana = Analysis(Ro, Bu, Lr, Ur, cyclonic=cyclonic)

    L = ana.L
    H = ana.H
    g = ana.g
    c02 = g*H
    k = ana.l
    f = ana.f
    omega = ana.omega

    x_axis = ana.x_axis
    t_axis = ana.time_axis
    v  = ana.u_wave
    u  = ana.v_wave
    eta  = ana.h_wave #no mean depth included
    N = 4
    T = 2*np.pi/omega
    Ttot = t_axis[-1]
    Tstart = Ttot - N*T  # when averaging begins
    istart = ((t_axis - Tstart)**2).argmin() 
    u2avg2D = u[istart:-1, :, :]
    v2avg2D = v[istart:-1, :, :]
    h2avg2D = eta[istart:-1, :, :]
    return u2avg2D, v2avg2D, h2avg2D, x_axis, c02, L

# u2avg2D, v2avg2D, h2avg2D, x_axis, c02, L = load_exp_vars(Ro, Bu, Lr, cyclonic)

def filter_wave(field):
    ''' expecting  dims (time, space, space)'''
    field_ft = fft.fftshift(fft.fft2(field, axes = (1, 2)), axes = (1, 2))
    nt, nx, ny = np.shape(field)
    field_ft[:, nx//2, :] = 0.
    filtered_field = fft.ifft2(fft.ifftshift(field_ft, axes = (1, 2)), axes = (1, 2))
    return np.real(filtered_field)



def linear_flux(u, v, eta, c02):
    factp2D = c02*eta
    Fxp2D, Fyp2D = np.mean(factp2D*u, 0), np.mean(factp2D*v, 0)
    return Fxp2D, Fyp2D

# Fxp2D, Fyp2D =  linear_flux(u2avg2D, v2avg2D, h2avg2D)

def interpF(Fxp2D, Fyp2D,x_axis, L ):
    IntrpFx = RGI((x_axis/L, x_axis/L), Fxp2D.T, method='linear')
    IntrpFy = RGI((x_axis/L, x_axis/L), Fyp2D.T, method='linear')
    return IntrpFx, IntrpFy

# IntrpFx, IntrpFy = interpF(Fxp2D, Fyp2D,x_axis, L )

# Dist = 3. # How far the integration contours are

def flux_arc(Fxp2D, Fyp2D,IntrpFx,IntrpFy, Dist):
     
    Angles = np.linspace(-np.pi/2, np.pi/2, 128)
    x_ArcCirc = Dist*np.cos(Angles)
    y_ArcCirc = Dist*np.sin(Angles)
    ArcCirc = [[x_ArcCirc[ii], y_ArcCirc[ii]] for ii in range(len(Angles))]
    caxsc = np.amax(abs(Fxp2D))
    FlxDensArc = np.cos(Angles)*IntrpFx(ArcCirc) + np.sin(Angles)*IntrpFy(ArcCirc)
    return  np.trapz(FlxDensArc, Dist*Angles)  # Not multiplying by L because we are only comparing

def flux_line(Fxp2D, Fyp2D, Dist, x_axis, L):
    iDist = ((x_axis/L + Dist)**2).argmin()
    FlxIn = np.trapz(Fxp2D[iDist:(1-iDist), iDist], x_axis[iDist:(1-iDist)]/L)
    return FlxIn


def flux_stack(Ro, Bu, Lr, cyclonic,  Ur =1000., filtered_wave = False):
    '''returns flux on larc cirlce and incoming side'''
    u2avg2D, v2avg2D, h2avg2D, x_axis, c02, L = load_exp_vars(Ro, Bu, Lr, cyclonic)
    if filtered_wave:
        u2avg2D = filter_wave(u2avg2D)
        v2avg2D = filter_wave(v2avg2D)
        h2avg2D = filter_wave(h2avg2D)

    Fxp2D, Fyp2D =  linear_flux(u2avg2D, v2avg2D, h2avg2D, c02)
    IntrpFx, IntrpFy = interpF(Fxp2D, Fyp2D, x_axis, L )
    Dist = 2. # How far the integration contours are
    iDist = ((x_axis/L + Dist)**2).argmin()# index where x=-2
    
    FlxArc = flux_arc(Fxp2D, Fyp2D,IntrpFx,IntrpFy, Dist)  
    FlxIn = flux_line(Fxp2D, Fyp2D, Dist, x_axis, L)
    print (FlxIn/FlxArc, FlxIn, FlxArc, filtered_wave)

    return FlxIn, FlxArc

def flux_ratio(Ro, Bu, Lr, cyclonic):
    FlxIn, FlxArc = flux_stack(Ro, Bu, Lr, cyclonic)
    FlxIn_scat, FlxArc_scat = flux_stack(Ro, Bu, Lr, cyclonic, filtered_wave =True)
    return FlxArc_scat/FlxIn

