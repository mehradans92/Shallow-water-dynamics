import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from ScientificCbar import ScientificCbar
from matplotlib.collections import LineCollection

matplotlib.rcParams.update({'font.size': 20})

## Loading data:
simulation_data = np.load('simulation_results.npz')

times = simulation_data['times']
x     = simulation_data['x']

H0 = simulation_data['H0']

g0 = simulation_data['g0']

simulated_h = simulation_data['simulated_h']
simulated_u = simulation_data['simulated_u']
simulated_u[:,:] = np.round(simulated_u[:,:],6)
Eta_B = simulation_data['Eta_B']

## Conservation of energy
dx = x[1] - x[0]
KE = np.sum(0.5 * simulated_h * simulated_u**2, axis=-1) * dx
PE = np.sum(0.5 * g0 * (simulated_h + Eta_B)**2, axis=-1) * dx
E_total = KE + PE

fig, axes = plt.subplots(2, 1, sharex=True, figsize=(10,6))


## Determine which times we want to show
Nt = len(times)

start = 0                    # start at the beginning
stop  = Nt//5                # only plot first fifth of simulation
skip  = (stop-start) // 10   # only use 10 lines to plot evolution

plot_inds = np.arange(start, stop, skip)

cmap = matplotlib.cm.viridis


## Plot h
h_lines = [np.column_stack([x, simulated_h[ii,:]-H0+Eta_B]) for ii in plot_inds]
h_segments = LineCollection(h_lines, cmap=cmap)
h_segments.set_array(times[plot_inds])
axes[0].add_collection(h_segments)
#axes[0].plot(x,Eta_B)



## Plot u
u_lines = [np.column_stack([x, simulated_u[ii,:]]) for ii in plot_inds]

u_segments = LineCollection(u_lines, cmap=cmap)
u_segments.set_array(times[plot_inds])
axes[1].add_collection(u_segments)

for ax in axes:
    ax.axis('tight')
    
axes[1].set_xlabel('x (m)')
axes[0].set_ylabel('$h -H0 +\eta_B$ (m)')
axes[1].set_ylabel('$u$ (m/s)')
    
# Add a colour bar showing time
cbar = plt.colorbar(h_segments, ax = axes)
ScientificCbar(cbar, label='time', units='(s)', labelpad=40, centre = False)

# Save to a png file
plt.savefig('line_plot.png', dpi=500)

fig, axes = plt.subplots(1, 2, sharey=True, figsize=(16,6))

# Plot the height field
cv  = np.max(np.abs(simulated_h - H0))
q0  = axes[0].pcolormesh(x, times, simulated_h-H0+Eta_B, vmin = -cv, vmax = +cv)
cb0 = plt.colorbar(q0, ax = axes[0])

# Plot the velocity field
cv  = np.max(np.abs(simulated_u))
q1  = axes[1].pcolormesh(x, times, simulated_u, vmin = -cv, vmax = +cv)
cb1 = plt.colorbar(q1, ax=axes[1])

# Labels
axes[0].set_xlabel('x (m)')
axes[1].set_xlabel('x (m)')
axes[0].set_ylabel('t (s)')

ScientificCbar(cb0, units = 'm',   label = 'Height Deviation $(h-H_0+\eta_B)$', labelpad = 40)
ScientificCbar(cb1, units = 'm/s', label = 'Velocity $(u)$',            labelpad = 40)

# Save to a png file
plt.savefig('hovmoller_plot.png', dpi=500)

fig, axes = plt.subplots(1, 1, sharey=True, figsize=(16,6))
plt.plot( times, KE - KE[0], label='KE' )
plt.plot( times, PE - PE[0], label='PE' )
plt.plot( times, E_total - E_total[0], label='KE + PE' )
plt.xlabel('time (s)')
plt.ylabel('Deviation from Initial Values')
plt.legend(loc='best')
plt.savefig('Cons_energy.png', dpi=500)