# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Reads a binary file and plots a cut in the x-y plane.
# Execute in ipython with '%run plot_bindata.py'

# Path to the binary data
N = 10
path = './N_{:d}/'.format(N)

# Spatial resolution
NX = 128
NY = 128
NZ = 64
shape = (NX,NY,NZ)

x = np.linspace(0, 2*np.pi, NX)

#%% plot del colormap
#
#plt.figure()
#for i in range(1,12):
#    field = np.fromfile(path+'th.{:04d}.out'.format(i),dtype=np.float32).reshape(shape,order='F')
#    plt.clf()
#    plt.title('i = {:}'.format(i))
#    plt.imshow(np.transpose(field[:,NY//2,:]), cmap = 'jet', origin = 'lower')#, vmin = , vmax = )
#    plt.xlabel('x')
#    plt.ylabel('z')
#    plt.colorbar()
#    plt.pause(1)
#plt.show()

#%% Plotemos el colormap para la ultima iteración

plt.figure()
field = np.fromfile(path+'th.0011.out',dtype=np.float32).reshape(shape,order='F')
plt.title('Temperatura potencial en el plano $XZ$', fontsize = 15)
plt.imshow(np.transpose(field[:,NY//2,:]), cmap = 'jet', origin = 'lower')#, vmin = , vmax = )
plt.xlabel('Grilla en x', fontsize = 15)
plt.ylabel('Grilla en z', fontsize = 15)
plt.colorbar()
plt.tight_layout()
plt.show()


#%% plot de temp(x)

plt.figure()
for i in range(1,12):
    field = np.fromfile(path+'th.{:04d}.out'.format(i),dtype=np.float32).reshape(shape,order='F')
    plt.clf()
    plt.title('i = {:}'.format(i))
    plt.plot(x, field[:, NY//2, NZ//2])
    plt.xlabel('x')
    plt.ylabel('temp')
    plt.grid(True)
    plt.pause(1)
plt.show()


#%% Ploteamos para la ultima iteracion

plt.figure()
field = np.fromfile(path+'th.0011.out'.format(i),dtype=np.float32).reshape(shape,order='F')
plt.title('')
plt.plot(x, field[:, NY//2, NZ//2], 'r', label = 'N = {:d}'.format(N))
plt.xlabel('Coordenada x', fontsize = 15)
plt.ylabel('Temperatura potencial', fontsize = 15)
plt.title('Corte de la temperatura potencial para $Y = NY/2$, $Z = NZ/2$', fontsize = 15)
plt.grid(True)
plt.legend()
plt.show()


