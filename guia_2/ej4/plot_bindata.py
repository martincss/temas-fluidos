# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Reads a binary file and plots a cut in the x-y plane.
# Execute in ipython with '%run plot_bindata.py'

# Path to the binary data
path = './N_10/'

# Spatial resolution
NX = 128
NY = 128
NZ = 64
shape = (NX,NY,NZ)

x = np.linspace(0, 2*np.pi, NX)

#%% plot del colormap
"""
plt.figure()
for i in range(1,12):
    field = np.fromfile(path+'th.{:04d}.out'.format(i),dtype=np.float32).reshape(shape,order='F')
    plt.clf()
    plt.title('i = {:}'.format(i))
    plt.imshow(np.transpose(field[:,NY//2,:]), cmap = 'jet', origin = 'lower')#, vmin = , vmax = )
    plt.xlabel('x')
    plt.ylabel('z')
    plt.colorbar()
    plt.pause(1)
plt.show()
"""
#%% plot de temp(x)

plt.figure()
for i in range(11,12):
    field = np.fromfile(path+'th.{:04d}.out'.format(i),dtype=np.float32).reshape(shape,order='F')
    plt.clf()
    plt.title('i = {:}'.format(i))
    plt.plot(x, field[:, NY//2, NZ//2])
    plt.xlabel('x')
    plt.ylabel('temp')
    plt.grid(True)
    plt.pause(1)
#plt.show()

#%% Fitteo
x_fit = x[:20]
temp_fit = field[:20, NY//2, NZ//2]
temp_fit /= np.mean(np.abs(temp_fit))
f = lambda x, A, lon, ph: A*np.sin(2*np.pi/lon*x+ph)

p_guess = [1.5, 0.8, -0.2]
popt, _ = curve_fit(f, x_fit, temp_fit, p0 = p_guess)
fit = f(x_fit, *popt)

plt.figure()
plt.plot(x_fit, fit, 'r')
plt.scatter(x_fit, temp_fit, color = 'b')
plt.show()

