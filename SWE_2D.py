import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

import math

# Domain parameters
Nx = 100
Ny = Nx # Assume unifrom spacing
Lx = 6.
Ly = 6.
dx = Lx / Nx
dy = Ly / Ny

# Creating the spatial computational grid:
x = np.arange(dx/2, Lx, dx) - Lx/2
y = np.arange(dy/2, Ly, dy) - Ly/2
# Physical parameters
f_c = 1 # coriolis force
H0 = 1.
g0 = 9.81
Eta_B = 0.2 * (np.exp( - (2*x / Lx)**2 ) * np.exp( - (2*y /Ly )**2 )) * np.ones((Nx,Ny))   # Bottom topography

# Gravity wave speed
c0 = np.sqrt(g0 * H0)

# Spatial Differentiation function (second-order)
def ddx(f):
    xp1 = np.roll(f, -1,axis=1)
    xm1 = np.roll(f,  1,axis=1)
    
    d = (xp1 - xm1) / (2*dx)
    
    return d

def ddy(f):
    
    yp1 = np.roll(f, -1, axis=0)
    ym1 = np.roll(f,  1, axis=0)
    
    g = (yp1 - ym1) / (2*dy)
    
    return g

# Spatial Differentiation function (forth-order)
def ddx4(f):
    
    xp1 = np.roll(f, -1,axis=1)
    xp2 = np.roll(f, -2,axis=1)
    xm1 = np.roll(f,  1,axis=1)
    xm2 = np.roll(f,  2,axis=1)
    
    dx4 = ( -1/12 * xp2 + 2/3 * xp1 - 2/3 * xm1 + 1/12 * xm2) / (dx)
    # These coefficients are based on the output of Coefficients.py
    return dx4

def ddy4(f):
    
    yp1 = np.roll(f, -1,axis=0)
    yp2 = np.roll(f, -2,axis=0)
    ym1 = np.roll(f,  1,axis=0)
    ym2 = np.roll(f,  2,axis=0)
    dy4 = ( -1/12 * yp2 + 2/3 * yp1 - 2/3 * ym1 + 1/12 * ym2) / (dy)
    # These coefficients are based on the output of Coefficients.py
    return dy4

    
### Constrain temporal grid:
# Initial conditions

# CASE a
u_init = np.zeros((Nx,Ny))
v_init = np.zeros((Nx,Ny))
X, Y = np.meshgrid(x, y)
def f(x,y):
    return  ( (1e-2* np.exp( - (x/ 0.2)**2 - (y / 0.2)**2))) 
h_init = H0 + f(X,Y) - Eta_B
ax = plt.axes(projection='3d')
ax.contour3D(X, Y, h_init - H0 + Eta_B, 100, cmap='viridis')
plt.savefig('h_initial.png', dpi=500)
# CASE b
# pi=math.pi
# u_init =   3e-2 * np.cos( 4*pi*x/Lx) + 1e-1 * np.exp( - (x / 0.2)**2 )
# h_init = H0 *np.ones(Nx)

# Time-stepping cfl factor
cfl = 0.1

# Target end-time for the simulation
final_time = np.maximum(5 * (Lx/2) / c0 , 5 * (Ly/2) / c0)

# Initial time
time0 = 0.

# State how frequently we want to store outputs
out_freq = final_time / 200
Nouts    = int(final_time / out_freq) + 1
out_ind  = 1
# The times for the stored outputs
times = np.arange(Nouts) * out_freq


# An array to store the spatio-temporal evolution
simulated_h = np.zeros((Nouts, Nx,Ny))
simulated_u = np.zeros((Nouts, Nx,Ny))
simulated_v = np.ones((Nouts, Nx,Ny))
# Store the initial conditions
simulated_u[0,:,:] = u_init
simulated_v[0,:,:] = v_init
simulated_h[0,:,:] = h_init

# Initialize some working variables
time = time0
u = u_init.copy()
v = v_init.copy()
h = h_init.copy() 
# Loop through time until we get to the target date
cnt = 0

# Discretization Schemes
# Temporal scheme order (Choose between 1, 2 or 3)
t_order = 3
# Spatial scheme order (Choose between 2 or 4)
s_order = 4

# Storing dudt and dhdt in an array
dudt = np.zeros((Nx,Ny,t_order))
dvdt = np.zeros((Nx,Ny,t_order))
dhdt = np.zeros((Nx,Ny,t_order))

while time < final_time:

    # Get dt
    dt = np.min([cfl * dx / (np.max(np.abs(u)) + c0),cfl * dy / (np.max(np.abs(v)) + c0)])

    if s_order == 2:
        dudt[:,:,0] = f_c*v - u * ddx(u) -v *ddy(u) - g0 * ddx(h+Eta_B) 
        dvdt[:,:,0] = -f_c*u - u * ddx(v) -v *ddy(v) - g0 * ddy(h+Eta_B)
        dhdt[:,:,0] = - ddx(u * h) - ddy(v * h)
    elif s_order == 4:
        dudt[:,:,0] = f_c*v - u * ddx4(u) -v *ddy4(u)  - g0 * ddx4(h+Eta_B) 
        dvdt[:,:,0] = -f_c*u - u * ddx4(v) -v *ddy4(v) - g0 * ddy4(h+Eta_B)
        dhdt[:,:,0] = - ddx4(u * h) - ddy4(v * h)

    # Update fields (1st-order temporal)
    if t_order == 1:
        # Compute RHS of evolution equations and save
        dt_now = dt
        u += dudt[:,:,0] * dt_now
        v += dvdt[:,:,0] * dt_now
        h += dhdt[:,:,0] * dt_now
    if t_order == 2:
        #Update fields (2nd-order temporal)
        if cnt == 0:
            dt_now = dt/10.
            u += dudt[:,:,0] * dt_now
            h += dhdt[:,:,0] * dt_now
            v += dvdt[:,:,0] * dt_now
        else:
            dt_now = dt
            u += ( 1.5 * dudt[:,:,0] - 0.5 * dudt[:,:,-1] ) * dt_now
            h += ( 1.5 * dhdt[:,:,0] - 0.5 * dhdt[:,:,-1])  * dt_now
            v += ( 1.5 * dvdt[:,:,0] - 0.5 * dvdt[:,:,-1])  * dt_now
    if t_order == 3:
        #Update fields (3rd-order temporal)
        if cnt == 0:
            dt_now = dt/10.
            u += dudt[:,:,0] * dt_now
            h += dhdt[:,:,0] * dt_now
            v += dvdt[:,:,0] * dt_now
        if cnt == 1:
            dt_now = dt/10.
            u += ( 1.5 * dudt[:,:,0] - 0.5 * dudt[:,:,-1] ) * dt_now
            h += ( 1.5 * dhdt[:,:,0] - 0.5 * dhdt[:,:,-1])  * dt_now
            v += ( 1.5 * dvdt[:,:,0] - 0.5 * dvdt[:,:,-1])  * dt_now
        else:
            dt_now = dt
            u += ( 5/3 * dudt[:,:,0] - 5/6 * dudt[:,:,-1] + 1/6 * dudt[:,:,-2] ) * dt_now
            h += ( 5/3 * dhdt[:,:,0] - 5/6 * dhdt[:,:,-1] + 1/6 * dhdt[:,:,-2])  * dt_now
            v += ( 5/3 * dvdt[:,:,0] - 5/6 * dvdt[:,:,-1] + 1/6 * dvdt[:,:,-2])  * dt_now
            
        
    # Update the history of dudt and dhdt values
    dudt = np.roll(dudt, -1, axis=2)
    dvdt = np.roll(dvdt, -1, axis=2)
    dhdt = np.roll(dhdt, -1, axis=2)

    # Update time 
    time += dt_now
    # u_vec[cnt] = u[cnt]
    # h_vec[cnt] = h[cnt]
    cnt += 1

    
    # Output, if appropraite
    if (time >= out_ind * out_freq):
        simulated_h[out_ind,:,:] = h 
        simulated_u[out_ind,:,:] = u
        simulated_v[out_ind,:,:] = v

        if (out_ind % 10 == 0):
            print('Stored output {0:d} of {1:d}.'.format(out_ind, Nouts-1))

        out_ind += 1

## Conservation of energy
dx = x[1] - x[0]
KE = np.sum(0.5 * simulated_h * (simulated_u**2 + simulated_v**2), axis=-1) * dx
PE = np.sum(0.5 * g0 * (simulated_h + Eta_B)**2, axis=-1) * dx

np.savez('simulation_results.npz',
         times = times, 
         x = x, 
         y = y,
         simulated_h = simulated_h, 
         simulated_u = simulated_u,
         simulated_v = simulated_v,
         Eta_B = Eta_B, 
         kinetic = KE,
         potential = PE,
         total = KE + PE,
         H0 = H0, 
         g0 = g0,
         Lx = Lx)

#h_anim = plt3D.h_anim(X, Y, simulated_h, Nouts*dt, "h")

