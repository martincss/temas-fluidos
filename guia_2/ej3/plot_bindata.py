import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.measurements import center_of_mass
# Reads a binary file and plots a cut in the x-y plane.
# Execute in ipython with '%run plot_bindata.py'



# Spatial resolution
NX = 128
NY = 128
NZ = 128
delta_t = 5e-3
delta_z = 2*np.pi/NZ
shape = (NX,NY,NZ)
iters = 20

centro_mas = np.zeros((iters, 3))
centro_menos = np.zeros((iters, 3))
vel_centro_mas = np.zeros((iters-1, 3))
vel_centro_menos = np.zeros((iters-1, 3))

omega = range(0,25,5)
vel_media_centro_mas = np.zeros(len(omega))
vel_media_centro_menos = np.zeros(len(omega))

dvx_dz = np.zeros(iters)
dvx_dz_media = np.zeros(len(omega))

for j in range(len(omega)):
    # Path to the binary data
    path = './omega_{:d}/'.format(omega[j])
    print path
    # plt.figure()
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

        # plt.clf()
        # plt.title('i = {:}'.format(i))
        # plt.imshow(np.transpose(h[:,NY//2,:]), cmap = 'jet', vmin = -15, vmax = 15, origin = 'lower')
        # plt.xlabel('x')
        # plt.ylabel('z')
        # plt.colorbar()
        # plt.pause(0.05)
    # plt.show()
    vel_centro_mas = (np.roll(centro_mas[:,2],-1,0) - np.roll(centro_mas[:,2],1,0)) / 2*delta_t
    vel_centro_menos = (np.roll(centro_menos[:,2],-1,0) - np.roll(centro_menos[:,2],1,0)) / 2*delta_t

    vel_media_centro_mas[j] = np.mean(vel_centro_mas[2:-2])
    vel_media_centro_menos[j] = np.mean(vel_centro_menos[2:-2])

    # calculo la derivada de vx (notacion: vx = u) con respecto a z
    dvx_dz = (np.roll(vx,-1,0) - np.roll(vx,1,0)) / 2*delta_z
    dvx_dz_media[j] = np.mean(abs(dvx_dz))

plt.clf()
plt.title(r'Velocidad de $h$ en funcion de $\Omega$')
plt.plot(omega, vel_media_centro_mas, 'x', omega, vel_media_centro_menos, 'x')
plt.show()

plt.clf()
plt.title('$\\dfrac{\\partial v_x}{\\partial z}$ de $h$ en funcion de $\\Omega$')
plt.plot(omega, dvx_dz_media, 'x')
plt.show()
