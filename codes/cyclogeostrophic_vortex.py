

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 10:11:50 2020

@author: jeff
"""


from dedalus import public as de
import matplotlib.pyplot as plt
import numpy as np
#from scipy.interpolate import interp2d
import h5py



import sys
import pickle

exp_dir =  sys.argv[1]


par_dict = pickle.load(open(exp_dir + 'IC/parameters.pkl', 'rb'))

dx = par_dict['dx']
nx = par_dict['nx_v']
ny = par_dict['ny_v']
Lx = par_dict['Lx_v']
Ly = par_dict['Ly_v']
f = par_dict['f']
g = par_dict['g']
mu = par_dict['mu']
l = par_dict['l']
H = par_dict['H']
L = par_dict['L']
hg = par_dict['hg']
cyclonic = par_dict['cyclonic']



xfbasis = de.Fourier('x', nx, interval=(-Lx/2,Lx/2), dealias=1) 
yfbasis = de.Fourier('y', ny, interval=(-Ly/2,Ly/2), dealias=1) 
domain = de.Domain([xfbasis, yfbasis], grid_dtype=np.float64)

kx = domain.elements(0)
ky = domain.elements(1)

x = domain.grid(0)
y = domain.grid(1)

h = domain.new_field(name='hg')
u = domain.new_field(name='ug')
v = domain.new_field(name='vg')

#slices = domain.dist.grid_layout.slices(scales=(1,1))

def gauss2d(x, y, sigx, sigy):
    return np.exp(-(x/sigx)**2/2-(y/sigy)**2/2)

if cyclonic:
    h['g'] = -hg*gauss2d(x, y, L/np.pi, L/np.pi)
else:
    h['g'] = hg*gauss2d(x, y, L/np.pi, L/np.pi)

#h['g'] = gauss2d(x, y, L/np.pi, L/np.pi)
#h['g'] =  hf.get('geoH')

h.differentiate('y', out=u)
u['g'] *= -1*g/f
h.differentiate('x', out=v)
v['g'] *= g/f

un = domain.new_field(name='un')
vn = domain.new_field(name='vn')
un1 = domain.new_field(name='un1')
vn1 = domain.new_field(name='vn1')

un['g'] = np.copy(u['g'])
vn['g'] = np.copy(v['g'])

un1['g'] = 0.
vn1['g'] = 0.

def error(un1, vn1, un, vn):
    resu = un1['g'] - un['g']
    resv = vn1['g'] - vn['g']
    
    res = np.max(np.sqrt(resu**2 +resv**2))
    return res

def uiter(un, vn):
    unx = domain.new_field(name='unx')
    uny = domain.new_field(name='uny')
    vnx = domain.new_field(name='vnx')
    vny = domain.new_field(name='vny')
    
    un.differentiate('x', out=unx)
    un.differentiate('y', out=uny)
    vn.differentiate('x', out=vnx)
    vn.differentiate('y', out=vny)
    
    uiter = u['g'] - (un['g']*vnx['g'] + vn['g']*vny['g'])/f
    viter = v['g'] + (un['g']*unx['g'] + vn['g']*uny['g'])/f
    return uiter, viter
    
    
resn1 = error(un1, vn1, un, vn)

print (resn1, 'res0')

resn = resn1 +1
i = 0

threshhold =  0.0001
while (resn1 > threshhold) and (resn1<resn):
    resn  = np.copy(resn1)
    un1['g'], vn1['g'] = uiter(un, vn)
    resn1 = error(un1, vn1, un, vn)
    
    if (resn1 > threshhold) and (resn1<resn):
        un['g'] = np.copy(un1['g'])
        vn['g'] = np.copy(vn1['g'])
    
    print(i, resn1, resn)
    i+=1



print (np.max(un['g']), np.max(u['g']))
extent = np.array([-Ly/2, Ly/2, -Lx/2, Lx/2])
plot = False
if plot:
    plt.subplot(2,2,1)
    plt.imshow(u['g'])
    plt.colorbar(orientation  = 'horizontal')
    plt.title('orig u')
    
    plt.subplot(2,2,2)
    plt.imshow(un['g'])
    plt.colorbar(orientation  = 'horizontal')
    plt.title('adjusted u')
    
    plt.subplot(2,2,3)
    plt.imshow(v['g'])
    plt.colorbar(orientation  = 'horizontal')
    plt.title('orig v')
    
    plt.subplot(2,2,4)
    plt.imshow(vn['g'])
    plt.colorbar(orientation  = 'horizontal')
    plt.title('adjusted v')
    plt.show()
    
    
    plt.subplot(2,1,1)
    plt.imshow(np.sqrt(u['g']**2 + v['g']**2))
    plt.colorbar(orientation  = 'horizontal')
    plt.title('orig')
    
    plt.subplot(2,1,2)
    plt.imshow(np.sqrt(un['g']**2 + vn['g']**2), extent = extent)
    plt.colorbar(orientation  = 'horizontal')
    plt.title('adjusted')
    plt.show()


#if cyclonic:
#    vortex
#else:
    

saveData=True
if saveData:   
    hf = h5py.File(exp_dir + 'VIC/cyclogeostrophic_vortex.h5', 'w')
    hf.create_dataset('geoH', data=h['g'])
    hf.create_dataset('geoU', data=un['g'])
    hf.create_dataset('geoV', data=vn['g'])
    hf.close()








