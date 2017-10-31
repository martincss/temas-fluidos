# -*- coding: utf-8 -*-
"""
Martín Carusso, Florencia Lazzari
Temas avanzados de dinámica de fluidos, 2c 2017
Guía 2

Programa para el análisis de los datos extraídos del código GHOST usando solver
ROTH (ejercicio 3)
"""

import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from scipy.ndimage.measurements import center_of_mass
# Reads a binary file and plots a cut in the x-y plane.
# Execute in ipython with '%run plot_bindata.py'



# Spatial resolution
NX = 128
NY = 128
NZ = 128
delta_t = 5e-3
delta_t_sampleo = 5e-2
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

derivada_media = np.zeros((len(omega), iters))

dvx_dz = np.zeros(iters)
dvx_dz_media = np.zeros(len(omega))

dvx_dy = np.zeros(iters)
dvx_dy_media = np.zeros(len(omega))

for j in range(len(omega)):
    # Path to the binary data
    path = './omega_{:d}/'.format(int(omega[j]))
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

        # para cada tiempo quiero un valor de derivada_vx_dz
        derivada_vx_dz = (np.roll(vx,-1,0) - np.roll(vx,1,0)) / 2*delta_z
        derivada_media[j,i] = np.mean(abs(derivada_vx_dz))

    vel_centro_mas = (np.roll(centro_mas[1:,2],-1,0) - np.roll(centro_mas[1:,2],1,0)) / 2*delta_t_sampleo
    vel_centro_menos = (np.roll(centro_menos[1:,2],-1,0) - np.roll(centro_menos[1:,2],1,0)) / 2*delta_t_sampleo

    # plt.clf()
    # plt.title('omega = {:d}'.format(int(omega[j])))
    # plt.plot(vel_centro_mas[1:-3], 'x')
    # plt.show()

    vel_media_centro_mas[j] = np.mean(vel_centro_mas[1:-3])
    vel_media_centro_menos[j] = np.mean(vel_centro_menos[1:-3])

    # calculo la derivada de vx (notacion: vx = u) con respecto a z
    dvx_dz = (np.roll(vx,-1,0) - np.roll(vx,1,0)) / 2*delta_z
    dvx_dz_media[j] = np.mean(abs(dvx_dz))

    dvx_dy = (np.roll(vx,-1,1) - np.roll(vx,1,1)) / 2*delta_y
    dvx_dy_media[j] = np.mean(abs(dvx_dy))


plt.clf()
plt.title(r'Velocidad de $h$ en funcion de $\Omega$', fontsize='15')
plt.plot(omega/100., vel_media_centro_mas, 'x', omega/100., vel_media_centro_menos, 'x')
plt.xlabel(r'$\Omega$', fontsize='15')
plt.ylabel(r'$V_h$', fontsize='15')
plt.xticks(fontsize='15')
plt.yticks(fontsize='15')
plt.show()

plt.clf()
plt.plot(omega/100., dvx_dz_media,'x', label=r'$\frac{\partial v_x}{\partial z}$')
plt.plot(omega/100., dvx_dy_media, 'x', label=r'$\frac{\partial v_x}{\partial y}$')
plt.legend()
plt.xlabel(r'$\Omega$', fontsize='15')
plt.xticks(fontsize='15')
plt.yticks(fontsize='15')
plt.show()

cociente = dvx_dz_media/dvx_dy_media

plt.clf()
plt.plot(omega/100., cociente, 'x')
plt.xlabel(r'$\Omega$', fontsize='15')
plt.ylabel(r'$\frac{\partial v_x}{\partial z} / \frac{\partial v_x}{\partial y}$', fontsize='15')
plt.xticks(fontsize='15')
plt.yticks(fontsize='15')
plt.show()

tiempo = np.arange(iters)*delta_t_sampleo

plt.clf()
plt.plot(tiempo[1:], derivada_media[0,1:], label=r'$\Omega$ = {:d}'.format(int(omega[0]/100)))
plt.plot(tiempo[1:], derivada_media[4,1:], label=r'$\Omega$ = {:d}'.format(int(omega[4]/100)))
plt.plot(tiempo[1:], derivada_media[8,1:], label=r'$\Omega$ = {:d}'.format(int(omega[8]/100)))
plt.plot(tiempo[1:], derivada_media[12,1:], label=r'$\Omega$ = {:d}'.format(int(omega[12]/100)))
plt.plot(tiempo[1:], derivada_media[20,1:], label=r'$\Omega$ = {:d}'.format(int(omega[20]/100)))
plt.xlabel('t', fontsize='15')
plt.ylabel(r'$<|\frac{\partial v_x}{\partial z}|>$', fontsize='15')
plt.xticks(fontsize='15')
plt.yticks(fontsize='15')
plt.legend()
plt.show()
