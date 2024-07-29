"""
Python script to create the setup figure
"""
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, Patch, Polygon
from AnalysisClass import Analysis

Ro = 0.03
Bu = 1.0
Lr = 4.0
Ur = 1000.

ana = Analysis(Ro, Bu, Lr, Ur, cyclonic=False)

# hv = ana.h[0]
# hw = ana.h[-1] - hv

plt.rcParams.update({  # LaTeX fonts
    "text.usetex": True,
    "font.family": "Helvetica",
    "font.size": 12,
    'text.latex.preamble': r"\usepackage{amsmath}"
})
# need the line below to use \boldsymbol (and other amsmath stuff)
#matplotlib.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]

nx = 128
L = 0.1
lmbda = 0.025
Uw = 0.5  # relative amplitude of the wave
xwf = -0.25  # position of the wave front
dwf = 0.01  # width of the wave front
xF, dF = -0.5, lmbda  # start (in x) and width of wave maker
xS, dS = 0.5 - lmbda, 1*dF  # start (in x) and width of sponge layer

k = 2*np.pi/lmbda
x = np.linspace(-0.5, 0.5, nx)
y = 1*x
X, Y = np.meshgrid(x, y)

# Fake data below. You should replace by actual model output
hv = np.exp(-0.5*(X**2 + Y**2)/(L/np.pi)**2)
hw = Uw*np.cos(k*X) * 0.5 * (1 + np.tanh(-(X - xwf)/dwf)) \
    * 0.5 * (1 + np.tanh((X + 0.5 - lmbda/2)/dwf))
# below, I zero out hv because I replaced vortex with streamplot
    
hvm = ana.h[0]
hw = ana.h[-1] - hvm
hw = hw/np.max(np.abs(hw))
vmin = np.min(hw)
vmax = np.max(hw)
# hw_n = hw[::4, ::4]
hw_n = hw
#
eta =  hw_n  # replace by actual data

alpha_op = 0.8

fig, ax = plt.subplots(1, 1)


# %% SSH
extent = [x[0], x[-1], y[0], y[-1]]
print (extent)
ax.imshow( eta, extent =extent, cmap='bwr', vmin=vmin, vmax=vmax, origin='lower')
#ax.contourf(X, Y, eta, 128, cmap='bwr', vmin=vmin, vmax=vmax)  # motion
ax.axhline(0., color='grey', linestyle=':')
ax.axvline(0., color='grey', linestyle=':')

# %% Domain
for rad in [ 0.03]:
    ax.add_patch(Circle((0.33, 0.32), radius=rad, clip_on=False, zorder=10,
                        linewidth=1.5, edgecolor='k',  facecolor='none'))

plt.scatter(0.33, 0.32, marker='x', color='k', s=100)

ax.text(0.3, 0.39, r'$f\boldsymbol{\hat z}$', color='k', bbox=dict(boxstyle="Round", fc="w", alpha=alpha_op, lw=0.), fontsize = 'x-large')

# %% Forcing and sponge layers
ax.add_patch(Rectangle((xF, -0.5), dF, 1., fill=False, hatch='--'))  # forcing
ax.add_patch(Rectangle((xS, -0.5), dS, 1., fill=False, hatch='x'))  # sponge

axTop = ax.secondary_xaxis('top')
# axTop.tick_params(axis='x', direction='inout')
axTop.set_xticks([xF+dF/2, xS+dS/2])
axTop.tick_params(width=0)
axTop.set_xticklabels(['Wave maker', 'Sponge layer'], fontsize = 'x-large')


# %% Wave
ax.annotate(r"$\boldsymbol{c_g} = \partial_k\omega\, \boldsymbol{\hat x}$",
            xy=(-0.35, 0.45), xytext=(-0.25, 0.4),
            bbox=dict(boxstyle="Round", fc="w", alpha=alpha_op, lw=0.),
            arrowprops=dict(arrowstyle='-'), fontsize = 'x-large')
for ya in [0.35, 0.25]:
    ax.annotate("", xy=(-0.35, ya), xytext=(-0.25, 0.4),
                arrowprops=dict(arrowstyle='-'))
# ax.text(-0.25, 0.4, r"$\mathbf{c_g} = \partial_k\omega\, \bf{\hat x}$",
#         bbox=dict(boxstyle="Round", fc="w", alpha=0.5, lw=0.))
for ya in [-0.45 + 0.1*ii for ii in range(10)]:
    ax.arrow(-0.42, ya, 0.05, 0., width=0.005)


# Wavelength markers

xls = [-0.3+3.75*lmbda, -0.3+4.75*lmbda]  # x-positions
yl = -.4
for xl in xls:
    ax.plot([xl]*2, [yl-.015, yl+.015], 'k')
ax.plot(xls, [yl]*2, 'k')
ax.annotate(r"$\lambda = 2\pi/k$", xy=(.5*(xls[0]+xls[1]), yl),
            xytext=(.5*(xls[0]+xls[1])+2*lmbda, yl+0.05),
            bbox=dict(boxstyle="Round", fc="w", alpha=alpha_op, lw=0.),
            arrowprops=dict(arrowstyle='-',
                            connectionstyle='angle3,angleA=0,angleB=90'),
                            fontsize  = 'x-large')


# # %% Vortex
# # u, v = Y*hv/L**2, -X*hv/L**2  # from geostrophic balance, all constants = 1
# u, v = Y*hv/(L/np.pi)**2, -X*hv/(L/np.pi)**2  
# #u = ana.v[0]/ np.max(np.abs(ana.v[0]))
# #v = ana.u[0]/ np.max(np.abs(ana.u[0]))
# # streamplot
# # seed_points = np.array([[0]*6, [-0.02-0.06*ii for ii in range(6)]])
# slx = slice(nx//4, 3*nx//4)
# #slx = slice(nx//2 -nx//8, nx//2 + nx//8)
# ax.streamplot(
#     X[slx, slx], Y[slx, slx], u[slx, slx], v[slx, slx],
#     color=(u[slx, slx]**2+v[slx, slx]**2)**.5, cmap='Greens',
#     # start_points=seed_points.T,
#     arrowstyle='fancy', density=0.5)

# for xl in [0, L]:
#     ax.plot([xl]*2, [0., 0.3], 'k--')
    
# ax.annotate("", xy=(0., 0.28), xytext=(L, 0.28),
#             arrowprops=dict(arrowstyle="<->"))
# ax.text(L/2, 0.29, "$L$", horizontalalignment='center',
#         verticalalignment='bottom')
# # ax.annotate("$L$", xy=(L/2, 0.3), xytext=(L/2, 0.35),
# #             bbox=dict(boxstyle="Round", fc="w", alpha=0.9, lw=0.),
# #             arrowprops=dict(arrowstyle='-'))

# for xl in [-L/2, L/2]:
#     ax.plot([xl]*2, [0., 0.3], 'k--')

# ax.annotate(r"$\lambda = 2\pi/k$", xy=(.5*(xls[0]+xls[1]), yl),
#             xytext=(.5*(xls[0]+xls[1])+2*lmbda, yl+0.05),
#             bbox=dict(boxstyle="Round", fc="w", alpha=0.5, lw=0.),
#             arrowprops=dict(arrowstyle='-',
#                             connectionstyle='angle3,angleA=0,angleB=90'))

# ax.annotate("", xy=(-L/2., 0.28), xytext=(L/2, 0.28),
#             arrowprops=dict(arrowstyle="<->"))
# ax.text(0, 0.18, "$L$", horizontalalignment='center',
#         verticalalignment='bottom', bbox=dict(boxstyle="Round", fc="w", alpha=alpha_op, lw=0.) , fontsize = 'x-large')
# ax.annotate("$L$", xy=(L/2, 0.3), xytext=(L/2, 0.35),
#             bbox=dict(boxstyle="Round", fc="w", alpha=0.9, lw=0.),
#             arrowprops=dict(arrowstyle='-'))



# %% sides labelling
plt.text(-0.57, 0.1, "Incoming side", rotation="vertical", fontsize = 'x-large')
# ax.set_ylabel('Incoming side')
axRight = ax.secondary_yaxis('right')
axRight.set_yticks([])
axRight.set_ylabel('Outgoing side', fontsize = 'x-large')


# %% axes
ax.set_aspect('equal', 'box')
# ax.tick_params(bottom=True, top=True, left=True)  # , right=False)
# ax.tick_params(labelbottom=True, labeltop=True,
#                labelleft=True)  # , labelright=False)
# ax.set_xticks([-0.5, 0., 0.5])
# ax.set_xticklabels(['$x=-L_x/2$', '$x=0$', '$x=L_x/2$'])
# ax.set_yticks([-0.5, 0., 0.5])
# ax.set_yticklabels(['$y=-L_x/2$', '$y=0$', '$y=L_x/2$'])

ax.set_xticks([-0.5, 0., 0.5])
ax.set_xticklabels(['$-L_x/2$', '$x=0$', '$L_x/2$'], fontsize = 'x-large')
ax.set_yticks([-0.5, 0., 0.5])
ax.set_yticklabels(['$-L_x/2$', '$y=0$', '$L_x/2$'], fontsize = 'x-large')



# plt.savefig('setup.eps')
# plt.savefig('setup.pdf')
# plt.savefig('setup.png', dpi=600)

zoom_range = 0.12207
line_size = 30
plt.plot([-zoom_range, zoom_range], [zoom_range, zoom_range], 'k-.', markersize = line_size)
plt.plot([-zoom_range, -zoom_range], [-zoom_range, zoom_range], 'k-.', markersize = line_size)
plt.plot([zoom_range, zoom_range], [zoom_range, -zoom_range], 'k-.', markersize = line_size)
plt.plot([-zoom_range, zoom_range], [-zoom_range, -zoom_range], 'k-.', markersize = line_size)

plt.show()
