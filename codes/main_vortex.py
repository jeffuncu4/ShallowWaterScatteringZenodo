from dedalus import public as de
import numpy as np
import h5py
import pickle
import sys

import matplotlib.pyplot as plt




exp_dir = sys.argv[1]

par_dict = pickle.load(open(exp_dir + 'IC/parameters.pkl', 'rb'))


nx = par_dict['nx_v']
ny = par_dict['ny_v']
Lx = par_dict['Lx_v']
Ly = par_dict['Ly_v']

#num_grid = 64
#dx = 500
#
#nx = num_grid
#ny = num_grid
#Lx = num_grid*dx
#Ly = num_grid*dx

f = par_dict['f']
g = par_dict['g']
mu = par_dict['mu']
l = par_dict['l']
H = par_dict['H']
tau_s = par_dict['tau_s']
num_iter_vortex = par_dict['num_iter_vortex']
dt = par_dict['dt']
save_iter_vortex = par_dict['save_iter_vortex']

xbasis = de.Fourier('x', nx, interval=(-Lx/2,Lx/2), dealias=3/2) 
ybasis = de.Fourier('y', ny, interval=(-Ly/2,Ly/2), dealias=3/2)

domain = de.Domain([xbasis, ybasis], grid_dtype=np.float64)
x = domain.grid(0)
y = domain.grid(1)


# Create problem and add constants
problem = de.IVP(domain, variables=['u', 'v', 'h'])
problem.parameters['g'] = g
problem.parameters['mu'] = mu
problem.parameters['f'] = f
problem.parameters['H'] = H
problem.parameters['l'] = l
problem.parameters['tau_s'] = tau_s*0.1
problem.parameters['pi'] = np.pi
problem.parameters['e'] = np.e
problem.parameters['T'] = np.pi*2/f


# import windows for forcing
slices = domain.dist.grid_layout.slices(scales=(1,1))

circularWindowData = h5py.File(exp_dir + 'VIC/circularWindow.h5', 'r')
circularWindow = domain.new_field(name = 'circular_window')
circularWindow['g'] = circularWindowData['circularWindow'][slices]




problem.parameters['cw'] = circularWindow


#problem.substitutions['dis(A)'] = "mu*(dx(dx(A)) + dy(dy(A)))" 
problem.substitutions['dis(A)'] = "-mu*d(A, x=4, y=4)" 
problem.substitutions['udg(A)'] = "u*dx(A) + v*dy(A)"

# Main equation, with linear terms on the LHS and nonlinear terms on the RHS
#problem.add_equation("dt(u) - f*v + g*dx(h) - dis(u) = -udg(u) \
#                        + (-u)/tau_s*cw*(1 - 1./(1 + e**(-20./T*(t -T/2))))")
#
#
#problem.add_equation("dt(v) + f*u + g*dy(h) - dis(v) = -udg(v) \
#                      + (-v)/tau_s*cw*(1 - 1./(1 + e**(-20./T*(t -T/2))))")
#
#
#problem.add_equation("dt(h) - dis(h) = -(h)*(dx(u) + dy(v)) \
#                     - (u)*dx(h) - (v)*dy(h) \
#                      + (H-h)/tau_s*cw*(1 - 1./(1 + e**(-20./T*(t -T/2))))")

problem.add_equation("dt(u) - f*v + g*dx(h) - dis(u) = -udg(u) \
                        + (-u)/tau_s*cw")


problem.add_equation("dt(v) + f*u + g*dy(h) - dis(v) = -udg(v) \
                      + (-v)/tau_s*cw")


problem.add_equation("dt(h) - dis(h) = -(h)*(dx(u) + dy(v)) \
                     - (u)*dx(h) - (v)*dy(h) \
                      + (H-h)/tau_s*cw")




solver = problem.build_solver(de.timesteppers.RK222)
u = solver.state['u']
v = solver.state['v']
h = solver.state['h']


# import geostorphic flow for initial conditions
geoData= h5py.File(exp_dir + 'VIC/cyclogeostrophic_vortex.h5', 'r')

h['g'] = H + geoData.get('geoH')[slices]
u['g'] = geoData.get('geoU')[slices]
v['g'] = geoData.get('geoV')[slices]


solver.stop_iteration = num_iter_vortex
solver.stop_sim_time = np.inf
solver.stop_wall_time = np.inf

analysis = solver.evaluator.add_file_handler(exp_dir + 'IC', iter=num_iter_vortex//100)
#analysis = solver.evaluator.add_file_handler(exp_dir + 'IC', iter=num_iter_vortex - 1)
analysis.add_system(solver.state, layout='g')

from dedalus.extras import flow_tools
 #CFL
CFL = flow_tools.CFL(solver, initial_dt=dt, cadence=1, safety=0.8, max_dt= dt,
                     max_change=1.5, threshold=0.05)
CFL.add_velocities(('u', 'v'))
CFLSwitch=True


import time

# Flow properties
flow_cadence    = 1
flow_property   = "h"
flow_name       = 'height'
flow            = flow_tools.GlobalFlowProperty(solver, cadence=flow_cadence)
flow.add_property(flow_property, name=flow_name)



# Main loop
start_time = time.time()

print (num_iter_vortex, 'num iter')
while solver.ok:
    if CFLSwitch:
        dt = CFL.compute_dt()
    solver.step(dt)
    
    if solver.iteration % 100 == 0:
        print('Completed iteration {}'.format(solver.iteration))
        if np.isnan(flow.max(flow_name)):
            raise ValueError('blew up')
#    if solver.iteration % 200 == 0:
#        plt.imshow(solver.state['h']['g'])
#        plt.colorbar()
#        plt.show()        

end_time = time.time()
print('Runtime:', end_time-start_time)




    


