import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.measurements import center_of_mass
# Reads a binary file and plots a cut in the x-y plane.
# Execute in ipython with '%run plot_bindata.py'

# Path to the binary data
path = './omega_10/'

# Spatial resolution
NX = 128
NY = 128
NZ = 128
delta_t = 5e-3
shape = (NX,NY,NZ)
iters = 21

centro_mas = np.zeros((iters, 3))
centro_menos = np.zeros((iters, 3))
vel_centro_mas = np.zeros((iters-1, 3))
vel_centro_menos = np.zeros((iters-1, 3))

plt.figure()
for i in range(1, iters):
    vx = np.fromfile(path+'vx.{:04d}.out'.format(i),dtype=np.float32).reshape(shape,order='F')
    vy = np.fromfile(path+'vy.{:04d}.out'.format(i),dtype=np.float32).reshape(shape,order='F')
    vz = np.fromfile(path+'vz.{:04d}.out'.format(i),dtype=np.float32).reshape(shape,order='F')

    wx = np.fromfile(path+'wx.{:04d}.out'.format(i),dtype=np.float32).reshape(shape,order='F')
    wy = np.fromfile(path+'wy.{:04d}.out'.format(i),dtype=np.float32).reshape(shape,order='F')
    wz = np.fromfile(path+'wz.{:04d}.out'.format(i),dtype=np.float32).reshape(shape,order='F')

    h = vx*wx + vy*wy + vz*wz
    h_mas = h*(h > 0)
    h_menos = h*(h < 0)
    centro_mas[i,0], centro_mas[i,1], centro_mas[i,2] = center_of_mass(h_mas)
    centro_menos[i,0], centro_menos[i,1], centro_menos[i,2] = center_of_mass(h_menos)


    plt.clf()
    plt.title('i = {:}'.format(i))
    plt.imshow(np.transpose(h[:,NY//2,:]), cmap = 'jet', vmin = -15, vmax = 15, origin = 'lower')
    plt.xlabel('x')
    plt.ylabel('z')
    plt.colorbar()
    plt.pause(0.05)
plt.show()


vel_centro_mas = (np.roll(centro_mas,-1,0) - np.roll(centro_mas,1,0)) / 2*delta_t
vel_centro_menos = (np.roll(centro_menos,-1,0) - np.roll(centro_menos,1,0)) / 2*delta_t
