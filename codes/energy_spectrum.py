from simulationClass import Simulation 
from AnalysisClass import Analysis
import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft


Ro = 0.04 # must make these floats for proper naming conventions
Bu = 1.0
Lr = 2.0
Ur = 1000.0
cyclonic  = True

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

    return u[istart], v[istart], eta[istart], x_axis, k

g = 9.81
def energy(hfft, ufft, vfft):
    return (hfft*ufft**2 + vfft**2 +g*hfft**2)/2.

u, v, eta, x_axis, k = load_exp_vars(Ro, Bu, Lr, cyclonic,  Ur =1000.)
print (np.shape(u))

ufft = fft.fft2(u)
vfft = fft.fft2(v)
etafft = fft.fft2(eta)

# Shift the zero frequency component to the center
ufft = fft.fftshift(ufft)
vfft = fft.fftshift(vfft)
etafft = fft.fftshift(etafft)


# Compute the wavenumber axes
dx = x_axis[1] - x_axis[0]
n = len(x_axis)
kx = fft.fftshift(fft.fftfreq(n, d=dx)) * 2 * np.pi/k
ky = kx  # Assuming square domain

# Create a meshgrid of the wavenumber axes
kx_grid, ky_grid = np.meshgrid(kx, ky)


# energy_spec_max = energy(etafft, ufft, vfft)[256, 286]
energy_spec_max = 2.8e4
mask = (ky_grid == 0)
ufft[mask] = 0
vfft[mask] = 0
etafft[mask] = 0

energy_spec = energy(etafft, ufft, vfft)
cmap = 'Blues'

radius = 1
theta = np.linspace(0, 2 * np.pi, 100)

# Parametric equations for the circle
xc = radius * np.cos(theta)
yc = radius * np.sin(theta)

label_size = 'x-large'

plt.figure(figsize = (5, 6))
plt.plot(xc, yc,'r--')
plt.imshow(np.abs(energy_spec), extent=[kx.min(), kx.max(), ky.min(), ky.max()],cmap =cmap)
# plt.imshow(np.abs(energy_spec),cmap =cmap)
cbar = plt.colorbar()
cbar.set_label('E')
plt.title('Scattered Wave Energy Spectrum')
plt.xlabel(r'$k/k_i$', fontsize = label_size)
plt.ylabel(r'$l/k_i$', fontsize = label_size)
plt.xlim(0,1.2)
plt.ylim(-1.2,1.2)
plt.show()