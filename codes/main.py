from dedalus import public as de
import numpy as np
import h5py
import pickle
import sys
import matplotlib.pyplot as plt


exp_dir =  sys.argv[1] 

par_dict = pickle.load(open(exp_dir + 'IC/parameters.pkl', 'rb'))


nx = par_dict['nx']
ny = par_dict['ny']
Lx = par_dict['Lx']
Ly = par_dict['Ly']
f = par_dict['f']
g = par_dict['g']
mu = par_dict['mu']
l = par_dict['l']
omega = par_dict['omega']
H = par_dict['H']
tau_s = par_dict['tau_s']
tau_w = par_dict['tau_w']
num_iter = par_dict['num_iter']
dt = par_dict['dt']
save_iter = par_dict['save_iter']
max_writes = par_dict['max_writes']
vortex_name = par_dict['vortex_name']


uw = par_dict['uw']
hw = par_dict['hw']
vw = par_dict['vw']


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


problem.parameters['u_force'] = uw
problem.parameters['v_force'] = vw
problem.parameters['h_force'] = hw
problem.parameters['l'] = l
problem.parameters['omega'] = omega
problem.parameters['tau_s'] = tau_s
problem.parameters['tau_w'] = tau_w
problem.parameters['pi'] = np.pi
problem.parameters['e'] = np.e
problem.parameters['T'] = np.pi*2/f


# import windows for forcing
slices = domain.dist.grid_layout.slices(scales=(1,1))
windowData = h5py.File(exp_dir + 'IC/forcingWindows.h5', 'r')

waveForcingWindow= domain.new_field(name = 'waveWindow')
waveForcingWindow.meta['x']['constant'] = True
waveForcingWindow['g'] = windowData.get('waveForcingWindow')[slices]

spongeWindow = domain.new_field(name = 'sponge_window')
spongeWindow.meta['x']['constant'] = True
spongeWindow['g'] = windowData['spongeWindow'][slices]



problem.parameters['sw'] = spongeWindow
problem.parameters['ww'] = waveForcingWindow



#problem.substitutions['dis(A)'] = "mu*(dx(dx(A)) + dy(dy(A)))" 
problem.substitutions['dis(A)'] = "-mu*d(A, x=4, y=4)" 
problem.substitutions['udg(A)'] = "u*dx(A) + v*dy(A)"

# Main equation, with linear terms on the LHS and nonlinear terms on the RHS. 

problem.add_equation("dt(u) - f*v + g*dx(h) - dis(u) = -udg(u) \
                      + (u_force*cos(l*y-omega*t) -u)/tau_w*ww + (-u)/tau_s*sw")

problem.add_equation("dt(v) + f*u + g*dy(h) - dis(v) = -udg(v) \
                     + (v_force*sin(l*y-omega*t) -v)/tau_w*ww + (-v)/tau_s*sw ")


problem.add_equation("dt(h) - dis(h) = -(h)*(dx(u) + dy(v)) \
                     - (u)*dx(h) - (v)*dy(h) \
                    + (H + h_force*sin(l*y-omega*t) -h)/tau_w*ww \
                     + (H-h)/tau_s*sw")
# Note that these equations use f>0, but x and y are swapped and this results
# in f<0 with the new definition of the axes.

solver = problem.build_solver(de.timesteppers.RK222)
u = solver.state['u']
v = solver.state['v']
h = solver.state['h']


# import geostorphic flow for initial conditions
geoData = h5py.File(exp_dir + 'IC/' + 'expanded_' + vortex_name , 'r')

h['g'] = geoData.get('geoH')[slices] 
u['g'] = geoData.get('geoU')[slices] 
v['g'] = geoData.get('geoV')[slices]

#h['g'] = H
#u['g'] = 0.
#v['g'] = 0.


solver.stop_iteration = num_iter
solver.stop_sim_time = np.inf
solver.stop_wall_time = np.inf


analysis = solver.evaluator.add_file_handler(exp_dir + 'data', iter=save_iter, max_writes=max_writes) # save by sim time
analysis.add_system(solver.state, layout='g')

uy = domain.new_field(name='uy')
vx = domain.new_field(name='vx')
v.differentiate('x', out=vx)
u.differentiate('y', out=uy)


curl = vx-uy


analysis.add_task(curl, layout='g', name='zeta')


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
print (num_iter, 'num iter')
while solver.ok:
    if CFLSwitch:
        dt = CFL.compute_dt()
    solver.step(dt)
    if solver.iteration % 100 == 0:
        print('Completed iteration {}/{}'.format(solver.iteration, num_iter))
        #print('Completed iteration {}'.format(solver.iteration) + ' ' + num_iter)
        if np.isnan(flow.max(flow_name)):
            raise ValueError('blew up')
#    if solver.iteration % 1000 == 0:
#        plt.imshow(solver.state['h']['g'])
#        plt.colorbar()
#        plt.show()

end_time = time.time()
print('Runtime:', end_time-start_time)




    


