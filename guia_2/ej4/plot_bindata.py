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
for i in range(1,2):
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
f = lambda x, A, lon, ph: A*np.sin(lon*x+ph)

popt, _ = curve_fit(f, x_fit, temp_fit)
fit = f(x_fit, *popt)

plt.plot(x_fit, fit, 'r')
plt.plot(x_fit, temp_fit, 'b')
plt.show()

