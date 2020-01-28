import numpy as np
import matplotlib as plt
import math

# Domain parameters
Nx = 512
Lx = 6.
dx = Lx / Nx

# Creating the spatial computational grid:
x = np.arange(dx/2, Lx, dx) - Lx/2

# Physical parameters
H0 = 1.
g0 = 9.81
Eta_B = 0.00005*(np.exp( - (2*x / Lx)**2 )) * np.ones(Nx)    # Bottom topography

# Gravity wave speed
c0 = np.sqrt(g0 * H0)

# Spatial Differentiation function (second-order)
def ddx(f):
    
    xp1 = np.roll(f, -1)
    xm1 = np.roll(f,  1)
    
    d = (xp1 - xm1) / (2*dx)
    
    return d

# Spatial Differentiation function (forth-order)
def ddx4(f):
    
    xp1 = np.roll(f, -1)
    xp2 = np.roll(f, -2)
    xm1 = np.roll(f,  1)
    xm2 = np.roll(f,  2)
    
    d4 = ( -1/12 * xp2 + 2/3 * xp1 - 2/3 * xm1 + 1/12 * xm2) / (dx)
    # These coefficients are based on the output of Coefficients.py
    return d4

    
### Constrain temporal grid:
# Initial conditions

# CASE a
u_init = np.zeros(Nx)
h_init = H0 +   1e-2 * np.exp( - (x / 0.2)**2 ) - Eta_B
#h_init = H0* np.ones(Nx) - Eta_B

# CASE b
# pi=math.pi
# u_init =   3e-2 * np.cos( 4*pi*x/Lx) + 1e-1 * np.exp( - (x / 0.2)**2 )
# h_init = H0 *np.ones(Nx) - Eta_B

# Time-stepping cfl factor
cfl = 0.2

# Target end-time for the simulation
final_time = 3 * (Lx/2) / c0

# Initial time
time0 = 0.

# State how frequently we want to store outputs
out_freq = final_time / 200
Nouts    = int(final_time / out_freq) + 1
out_ind  = 1
# The times for the stored outputs
times = np.arange(Nouts) * out_freq

# An array to store the spatio-temporal evolution
simulated_h = np.zeros((Nouts, Nx))
simulated_u = np.zeros((Nouts, Nx))
# Store the initial conditions
simulated_u[0,:] = u_init
simulated_h[0,:] = h_init

# Initialize some working variables
time = time0
u = u_init.copy()
h = h_init.copy() 
# Loop through time until we get to the target date
cnt = 0

# Discretization Schemes
# Temporal scheme order (Choose between 1, 2 or 3)
t_order = 3
# Spatial scheme order (Choose between 2 or 4)
s_order = 4

# Storing dudt and dhdt in an array
dudt = np.zeros((Nx,t_order))
dhdt = np.zeros((Nx,t_order))

while time < final_time:

    # Get dt
    dt = cfl * dx / (np.max(np.abs(u)) + c0)

    if s_order == 2:
        dudt[:,0] = - u * ddx(u) - g0 * ddx(h+Eta_B) 
        dhdt[:,0] = - ddx(u * h)
    elif s_order == 4:
        dudt[:,0] = - u * ddx4(u) - g0 * ddx4(h+Eta_B) 
        dhdt[:,0] = - ddx4(u * h)

    # Update fields (1st-order temporal)
    if t_order == 1:
        # Compute RHS of evolution equations and save
        dt_now = dt
        u += dudt[:,0] * dt_now
        h += dhdt[:,0] * dt_now
    if t_order == 2:
        #Update fields (2nd-order temporal)
        if cnt == 0:
            dt_now = dt/10.
            u += dudt[:,0] * dt_now
            h += dhdt[:,0] * dt_now
        else:
            dt_now = dt
            u += ( 1.5 * dudt[:,0] - 0.5 * dudt[:,-1] ) * dt_now
            h += ( 1.5 * dhdt[:,0] - 0.5 * dhdt[:,-1])  * dt_now
    if t_order == 3:
        #Update fields (3rd-order temporal)
        if cnt == 0:
            dt_now = dt/10.
            u += dudt[:,0] * dt_now
            h += dhdt[:,0] * dt_now
        if cnt == 1:
            dt_now = dt/10.
            u += ( 1.5 * dudt[:,0] - 0.5 * dudt[:,-1] ) * dt_now
            h += ( 1.5 * dhdt[:,0] - 0.5 * dhdt[:,-1])  * dt_now
        else:
            dt_now = dt
            u += ( 5/3 * dudt[:,0] - 5/6 * dudt[:,-1] + 1/6 * dudt[:,-2] ) * dt_now
            h += ( 5/3 * dhdt[:,0] - 5/6 * dhdt[:,-1] + 1/6 * dhdt[:,-2])  * dt_now
            
        
    # Update the history of dudt and dhdt values
    dudt = np.roll(dudt, -1, axis=1)
    dhdt = np.roll(dhdt, -1, axis=1)

    # Update time 
    time += dt_now
    # u_vec[cnt] = u[cnt]
    # h_vec[cnt] = h[cnt]
    cnt += 1

    
    # Output, if appropraite
    if (time >= out_ind * out_freq):
        simulated_h[out_ind,:] = h 
        simulated_u[out_ind,:] = u

        if (out_ind % 10 == 0):
            print('Stored output {0:d} of {1:d}.'.format(out_ind, Nouts-1))

        out_ind += 1



np.savez('simulation_results.npz', 
         times = times, 
         x = x, 
         simulated_h = simulated_h, 
         simulated_u = simulated_u,
         Eta_B = Eta_B, 
         H0 = H0, 
         g0 = g0,
         Lx = Lx)