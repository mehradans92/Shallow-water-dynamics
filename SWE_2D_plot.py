import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from matplotlib import gridspec
from zipfile import ZipFile
import os

mpl.rcParams.update({'font.size': 16})
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['xtick.major.pad']= 10

## Loading data:
simulation_data = np.load('simulation_results.npz')

times = simulation_data['times']
x     = simulation_data['x']
y     = simulation_data['y']
dx = x[1] - x[0]
dy = y[1] - y[0]
H0 = simulation_data['H0']
g0 = simulation_data['g0']

simulated_h = simulation_data['simulated_h']
simulated_u = simulation_data['simulated_u']
simulated_v = simulation_data['simulated_v']
Eta_B = simulation_data['Eta_B']

#KE = simulation_data['kinetic']
#PE = simulation_data['potential']

KE = np.sum(np.sum(0.5 * simulated_h * (simulated_u**2 + simulated_v**2), axis=-1) * dx,axis=-1)*dy
PE = np.sum(np.sum(0.5 * g0 * (simulated_h + Eta_B)**2, axis=-1) * dx, axis=-1)*dy
E_total = KE + PE
fig, axes = plt.subplots(1, 1, sharey=True, figsize=(16,6))
plt.plot( times, KE - KE[0], label='KE' )
plt.plot( times, PE - PE[0], label='PE' )
plt.plot( times, E_total - E_total[0], label='KE + PE' )
plt.xlabel('time (s)')
plt.ylabel('Deviation from Initial Values')
plt.legend(loc='best')
plt.savefig('Cons_energy.png', dpi=500)
X, Y = np.meshgrid(x, y)

#### PLOT Surface elevation only

# for i in range (0,len(simulated_h)):
#     Eta = simulated_h[i]-H0+Eta_B
#     Eta_0 = np.max(simulated_h[0]-H0+Eta_B)
#     #ax.contour3D(X, Y, simulated_h[i]-H0+Eta_B, 100, cmap='viridis')
#     #ax.plot_surface(X, Y, simulated_h[i]-H0+Eta_B, cmap = plt.cm.RdBu_r)
#     fig = plt.figure(figsize=(16,8))
#     ax = fig.gca(projection='3d')
#     surf = ax.plot_surface(X, Y, Eta, rstride = 1, cstride = 1,
#         cmap = plt.cm.RdBu_r, linewidth = 0, antialiased = True)
#     ax.set_zlim(-Eta_0,Eta_0)
#     ax.zaxis.set_rotate_label(False)
#     ax.view_init(30, 45)
#     ax.set_xlabel("x [m]", fontname = "serif", fontsize = 20,labelpad=18)
#     ax.set_ylabel("y [m]", fontname = "serif", fontsize = 20,labelpad=18)
#     ax.set_zlabel("Surface elevation [m]", fontname = "serif", fontsize = 20, labelpad=25, rotation=90)
#     ax.set_title(" $\eta(x,y,t)$ at $t = {:.2f} $ seconds".format(
#             times[i]), fontname = "serif", fontsize = 25, y=1.06)
#     #cbar = fig.colorbar(surf, shrink=0.5, aspect=6 )
#     #plt.ticks( fontsize = 16)
#     plt.savefig("Surface"+ str(i) + ".png", dpi=500)
#     plt.close()
#     if (i % 10 == 0):
#             print('Stored PNG snapshot {0:d} of {1:d}.'.format(i+1, len(simulated_h)-1))
# os.system('ffmpeg -r 20 -i Surface%01d.png -vcodec mpeg4 -y anim.mp4')
# os.system('rm Surface*.png')
# print ('Post-processing done!')


for i in range (0,len(simulated_h)):
    Eta_0 = np.max(simulated_h[0]-H0+Eta_B)
    fig, ax= plt.subplots(1, 2,figsize=(20,10))
    gs = gridspec.GridSpec(1, 2, width_ratios=[2.5, 1]) 
    [axi.set_axis_off() for axi in ax.ravel()]
    ## Plotting surface elevation profiles
    Eta = simulated_h[i]-H0+Eta_B
    ax0 = plt.subplot(gs[0],projection='3d')
    ax0.plot_surface(X, Y, Eta, rstride = 1, cstride = 1,
        cmap = plt.cm.RdBu_r, linewidth = 0, antialiased = True)
    ax0.set_zlim(-Eta_0,Eta_0)
    ax0.zaxis.set_rotate_label(False)
    ax0.view_init(30, 45)
    ax0.set_xlabel("x [m]", fontname = "serif", fontsize = 20,labelpad=18)
    ax0.set_ylabel("y [m]", fontname = "serif", fontsize = 20,labelpad=18)
    ax0.set_zlabel("Surface elevation [m]", fontname = "serif", fontsize = 20, labelpad=25, rotation=90)
    ax0.set_title(" $\eta(x,y,t)$ at $t = {:.2f} $ seconds".format(
            times[i]), fontname = "serif", fontsize = 25, y=1.06)

    ## Velocity profiles
    #ax1 = fig.add_subplot(1, 2, 2)
    ax1 = plt.subplot(gs[1])
    ax1.set_aspect(1)
    ax1.quiver(x, y, simulated_u[i],simulated_v[i],scale=0.1)
    ax1.set_title( "Velocity profile at $t = {:.2f} $ seconds".format(times[i]), fontname = "serif", fontsize = 25 , y=1.06)
    ax1.set_xlabel("x [m]", fontname = "serif", fontsize = 20)
    ax1.set_ylabel("y [m]", fontname = "serif", fontsize = 20)
    if (i % 10 == 0):
        print ('Stored PNG snapshot {0:d} of {1:d}.'.format(i+1, len(simulated_h)-1))
    plt.savefig("post"+ str(i) + ".png", dpi=200,bbox_inches='tight')
    plt.close()
    
os.system('ffmpeg -r 18 -i post%01d.png -vcodec mpeg4 -y anim.mp4')
#os.remove('.\post*.png')
print ('Post-processing done!') 
