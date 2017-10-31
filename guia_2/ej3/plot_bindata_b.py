import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.measurements import center_of_mass
# Reads a binary file and plots a cut in the x-y plane.
# Execute in ipython with '%run plot_bindata.py'



# Spatial resolution
NX = 128
NY = 128
NZ = 128
delta_t_sampleo = 5e-2
delta_t = 5e-3
delta_z = 2*np.pi/NZ
delta_y = 2*np.pi/NY
shape = (NX,NY,NZ)
iters = 20

centro_mas = np.zeros((iters, 3))
centro_menos = np.zeros((iters, 3))
vel_centro_mas = np.zeros((iters-1, 3))
vel_centro_menos = np.zeros((iters-1, 3))

omega = np.linspace(0,2000,17)
vel_media_centro_mas = np.zeros(len(omega))
vel_media_centro_menos = np.zeros(len(omega))

dvx_dz = np.zeros(iters)
dvx_dz_media = np.zeros(len(omega))

dvx_dy = np.zeros(iters)
dvx_dy_media = np.zeros(len(omega))

path = './omega_1000/'
for i in range(1, iters):
    vx = np.fromfile(path+'vx.{:04d}.out'.format(i),dtype=np.float32).reshape(shape,order='F')
    vy = np.fromfile(path+'vy.{:04d}.out'.format(i),dtype=np.float32).reshape(shape,order='F')
    vz = np.fromfile(path+'vz.{:04d}.out'.format(i),dtype=np.float32).reshape(shape,order='F')

    wx = np.fromfile(path+'wx.{:04d}.out'.format(i),dtype=np.float32).reshape(shape,order='F')
    wy = np.fromfile(path+'wy.{:04d}.out'.format(i),dtype=np.float32).reshape(shape,order='F')
    wz = np.fromfile(path+'wz.{:04d}.out'.format(i),dtype=np.float32).reshape(shape,order='F')

    if (i==1 or i==6 or i==12 or i==19):
        plt.clf()
        plt.title('t = {:}'.format(i*delta_t_sampleo), fontsize='20')
        plt.imshow(np.transpose(h[:,NY//2,:]), cmap = 'jet', vmin = -15, vmax = 15, origin = 'lower')
        plt.xlabel('x', fontsize='20')
        plt.ylabel('z', fontsize='20')
        plt.xticks(fontsize='20')
        plt.yticks(fontsize='20')
        cb = plt.colorbar()
        cb.ax.tick_params(labelsize=20)
        plt.pause(0.01)
        plt.show()
