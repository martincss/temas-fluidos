# -*- coding: utf-8 -*-
"""
Martín Carusso, Florencia Lazzari
Temas avanzados de dinámica de fluidos, 2c 2017
Guía 5

Programa para graficar la energía cinética, la energía potencial y la razón
entre ambas
(ejercicio 2)
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 17})
from scipy.interpolate import interp1d

###############################################################
# b)

#%% Importamos los datos de balance

#path = 'data\\'
path = 'data/'
file_balance = 'balance.txt'
tiempo, energia, enstrofia, epsilon = np.loadtxt(path + file_balance, unpack = True)
#energia = 0.5*energia
#enstrofia = 0.5*enstrofia

iters = 17
shape = 192*192*48
potencial = np.zeros(iters)

for t in range(0, iters):
    themp = np.fromfile(path+'th.{:04d}.out'.format(t+1),dtype=np.float32).reshape(shape,order='F')
    potencial[t] = np.mean(themp**2)

dt = 5e-3
sampleo_binario = 300
tiempo_potencial = np.arange(1, dt*sampleo_binario*iters, dt*sampleo_binario)

f = interp1d(tiempo_potencial, potencial)
potencial_inter = f(tiempo[2:-1])

#%% Graficamos estos datos

plt.figure()
plt.title(u'Energía')
plt.subplot(3,1,1)
plt.plot(tiempo[2:-1], energia[2:-1], color = 'r', label = u'Cinética $\\langle v^2 \\rangle$')
#plt.xlabel('Tiempo')
plt.xticks([0, 5, 10, 15, 20, 25], ()) # para que no me muestre el label en los
plt.ylabel(u'Energía Cinética')                  # de arriba
plt.grid(True)
plt.legend()

plt.subplot(3,1,2)
plt.plot(tiempo_potencial, potencial, color = 'b', label = u'Potencial $\\langle \\theta^2 \\rangle$')
#plt.xlabel('Tiempo')
plt.xticks([0, 5, 10, 15, 20, 25], ())
plt.ylabel(u'Energía Potencial')
plt.grid(True)
plt.legend()

plt.subplot(3,1,3)
plt.plot(tiempo[2:-1], energia[2:-1]/potencial_inter, color = 'g', label = u'Cociente $\\langle v^2 \\rangle/\\langle \\theta^2 \\rangle$')
plt.xlabel('Tiempo')
plt.ylabel('Cociente')
plt.grid(True)
plt.legend()


####################################################################
# c)
path = 'data/'

t_critico = np.where(tiempo == 15)[0][0]
suma_spec = np.zeros(97)
suma_para = np.zeros(25)
suma_perp = np.zeros(97)
i = 0

for t in range(len(tiempo[t_critico:])):
    i += 1

    # especto isotropo
    file_kspectrum = 'kspectrum.{:04d}.txt'.format(t+1+t_critico)
    k_spec = np.loadtxt(path + file_kspectrum, delimiter = '  ', usecols = (0,))
    spec = np.loadtxt(path + file_kspectrum, delimiter = '  ', usecols = (1,))
    suma_spec += spec

    # especto paralelo
    file_kspecpara = 'kspecpara.{:04d}.txt'.format(t+1+t_critico)
    k_para = np.loadtxt(path + file_kspecpara, delimiter = '  ', usecols = (0,))
    para = np.loadtxt(path + file_kspecpara, delimiter = '  ', usecols = (1,))
    suma_para += para

    # especto perpendicular
    file_kspecperp = 'kspecperp.{:04d}.txt'.format(t+1+t_critico)
    k_perp = np.loadtxt(path + file_kspecperp, delimiter = '  ', usecols = (0,))
    perp = np.loadtxt(path + file_kspecperp, delimiter = '  ', usecols = (1,))
    suma_perp += perp

mean_spec = suma_spec/i
mean_para = suma_para/i
mean_perp = suma_perp/i

#%% Graficamos todo

plt.figure()
plt.loglog(k_spec, mean_spec, color = 'b', label= u'$E(k)$')
plt.loglog(k_para, mean_para, color = 'g', label= '$E(k_{\\parallel})$')
plt.loglog(k_perp, mean_perp, color = 'r', label= '$E(k_{\\perp})$')
plt.xlabel('$k$', fontsize='20')
plt.ylabel('$E$', fontsize='20')
plt.title(u'Espectro de Energía')
plt.grid(True)
plt.legend()
plt.show()

####################################################################
# e)

t = np.where(tiempo_potencial == 21)[0][0]

Nx = 192
Ny = 192
Nz = 48
shape = (Nx,Ny,Nz)

Ny = 192
themperature = np.fromfile(path+'th.{:04d}.out'.format(t),dtype=np.float32).reshape(shape,order='F')
mean_themp = np.mean(themperature)
fluctuaciones = themperature - mean_themp

plt.clf()
plt.title('Temperatura a tiempo t = 21', fontsize='20')
plt.imshow(np.transpose(themperature[:,Ny//2,:]), cmap = 'jet', vmin = -0.6, vmax = 0.6, origin = 'lower')
plt.xlabel('x', fontsize='20')
plt.ylabel('z', fontsize='20')
plt.xticks(fontsize='20')
plt.yticks(fontsize='20')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=20)
plt.pause(0.01)
plt.show()

plt.clf()
plt.title('Fluctuaciones de Temperatura a tiempo t = 23.5', fontsize='20')
plt.imshow(np.transpose(fluctuaciones[:,Ny//2,:]), cmap = 'hot', vmin = -0.6, vmax = 0.6, origin = 'lower')
plt.xlabel('x', fontsize='20')
plt.ylabel('z', fontsize='20')
plt.xticks(fontsize='20')
plt.yticks(fontsize='20')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=20)
plt.pause(0.01)
plt.show()
