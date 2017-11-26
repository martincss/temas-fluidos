# -*- coding: utf-8 -*-
"""
Martín Carusso, Florencia Lazzari
Temas avanzados de dinámica de fluidos, 2c 2017
Guía 3

Programa para el análisis de los datos extraídos del código GHOST usando solver
BOUSS (guia 3, ejercicio 1)
"""
from __future__ import division

import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
# Reads a binary file and plots a cut in the x-y plane.
# Execute in ipython with '%run plot_bindata.py'

# Resolution
NX = 256
NY = 16
NZ = 128
shape = (NX,NY,NZ)
delta_t_sampleo = 0.5
iters = 21

delta_t_sampleo = 0.5
delta_z = 2*np.pi/NZ

z = np.arange(0,2*np.pi, delta_z)

# Estudio la evolucion temporal de u, ωy y de la temperatura

cant_N = 7

N = np.arange(1,cant_N)
rich_min = np.zeros(len(N))
vel_inest = np.zeros(len(N))

derivada_media = np.zeros((iters, NZ))
wy_media = np.zeros((iters, NZ))

index_pos_l =  np.zeros(iters)
index_pos_r = np.zeros(iters)
index_neg_l =  np.zeros(iters)
index_neg_r = np.zeros(iters)

frente_pos_l = np.zeros(iters)
frente_pos_r = np.zeros(iters)
frente_neg_l = np.zeros(iters)
frente_neg_r = np.zeros(iters)

pendiente_inest_media = 0
r_min = lambda N, u_0, gamma: N**2 / (u_0*gamma)**2

gamma = 10
u_0 = 1

def find(a, value):
    dif = np.abs(a-value)
    indice = np.argmin(dif)
    dif[indice] = np.max(dif)
    indice2 = np.argmin(dif)
    return np.min([indice, indice2]), np.max([indice, indice2])


for n in range(len(N)):
    # Path to the binary data
    path = './N{:d}_gamma10/'.format(int(N[n]))
    print path
    for t in range(0, iters):
        vx = np.fromfile(path+'vx.{:04d}.out'.format(t+1),dtype=np.float32).reshape(shape,order='F')
        wy = np.fromfile(path+'wy.{:04d}.out'.format(t+1),dtype=np.float32).reshape(shape,order='F')
        themp = np.fromfile(path+'th.{:04d}.out'.format(t+1),dtype=np.float32).reshape(shape,order='F')

        # plt.clf()
        # plt.title('t = {:}'.format((t)*delta_t_sampleo), fontsize='20')
        # plt.imshow(np.transpose(vx[:,NY//2,:]), cmap = 'jet', origin = 'lower')
        # plt.xlabel('x', fontsize='20')
        # plt.ylabel('z', fontsize='20')
        # plt.xticks(fontsize='20')
        # plt.yticks(fontsize='20')
        # cb = plt.colorbar()
        # cb.ax.tick_params(labelsize=20)
        # plt.pause(0.01)

        derivada_dvx_dz = (np.roll(vx,-1,2) - np.roll(vx,1,2)) / 2*delta_z
        derivada_media[t,:] = np.mean(derivada_dvx_dz, axis=(0,1))    # hay una derivada_media por instante de tiempo, el valor medio es espacial

        derivada_media_pos = derivada_media*(derivada_media>0)
        derivada_media_neg = derivada_media*(derivada_media<0)

        index_pos_l[t], index_pos_r[t] = find(derivada_media_pos[t,:], np.max(derivada_media[t, :])/4)
        index_neg_l[t], index_neg_r[t] = find(derivada_media_neg[t,:], np.min(derivada_media[t, :])/4)

        # plt.clf()
        # plt.title('t = {:}'.format(t*delta_t_sampleo), fontsize='20')
        # plt.plot(z, derivada_media[t,:])
        # plt.scatter(z[index_pos_l[t]], derivada_media[t,index_pos_l[t]], color='r', s=10)
        # plt.scatter(z[index_pos_r[t]], derivada_media[t,index_pos_r[t]], color='r', s=10)
        # plt.scatter(z[index_neg_l[t]], derivada_media[t,index_neg_l[t]], color='r', s=10)
        # plt.scatter(z[index_neg_r[t]], derivada_media[t,index_neg_r[t]], color='r', s=20)
        # plt.xlim([0,2*np.pi])
        # plt.ylim([-0.03,0.03])
        # plt.xlabel('z', fontsize='30')
        # plt.ylabel(r'$\partial v_x/\partial z$', fontsize='30')
        # plt.grid(True)
        # plt.xticks(fontsize='20')
        # plt.yticks(fontsize='20')
        # plt.pause(0.01)
        # plt.show()


        frente_pos_l[t] += z[index_pos_l[t]]
        frente_pos_r[t] += z[index_pos_r[t]]
        frente_neg_l[t] += z[index_neg_l[t]]
        frente_neg_r[t] += z[index_neg_r[t]]

    # plt.show()

    tiempo = np.arange(iters)*delta_t_sampleo

    # # esto no funciona porque la funcion es escalonada
    # derivada_index = np.diff(pos_frente) / delta_t_sampleo

    # # grafico uno de los frentes para que se vea esto
    # plt.clf()
    # plt.plot(tiempo, frente_pos_l, 'x')
    # plt.show()

    # en lugar de calcular la derivada a mano la voy a calcular a partir de
    # un ajuste lineal

    grupo_frentes=[frente_pos_l, frente_pos_r, frente_neg_l, frente_neg_r]

    pendiente_inest = np.zeros(len(grupo_frentes))
    i = 0
    for frente in grupo_frentes:

        params = np.polyfit(tiempo, frente, 1)

        y = lambda x, p: p[0] * x + p[1]
        fit = y(tiempo, params)

        # plt.clf()
        # plt.plot(tiempo, frente, 'x')
        # plt.plot(tiempo, fit)
        # plt.pause(0.01)

        pendiente_inest[i] = abs(params[0])
        i += 1

    # plt.show()

    vel_inest[n] = np.mean(pendiente_inest)
    rich_min[n] = r_min(N[n], u_0, gamma)

# plt.clf()
# plt.title('Velocidad de propagacion de la inestabilidad', fontsize='25')
# plt.plot(rich_min, vel_inest, 'ob')
# plt.xlabel(r'$Ri_{min}$', fontsize='25')
# plt.ylabel(r'$velocidad$', fontsize='25')
# plt.grid(True)
# plt.xticks(fontsize='25')
# plt.yticks(fontsize='25')
# plt.show()



gamma_array = np.arange(5,21,2)
long_onda_max = np.zeros(len(gamma_array))
long_onda_min = np.zeros(len(gamma_array))

from scipy.signal import argrelextrema

max_locales = []
min_locales = []

for g in range(len(gamma_array)):
    # Path to the binary data
    path = './N2_gamma{:d}/'.format(int(gamma_array[g]))
    print path
    # estudio unicamente para el ultimo instante de tiempo
    t = iters-1
    wy = np.fromfile(path+'wy.{:04d}.out'.format(t+1),dtype=np.float32).reshape(shape,order='F')
    wy_media[t,:] = np.mean(wy, axis=(0,1))

    wy_media_pos = wy_media[t,:]*(wy_media[t,:]>0)
    wy_media_neg = wy_media[t,:]*(wy_media[t,:]<0)

    x = argrelextrema(wy_media_pos, np.greater)
    y = argrelextrema(wy_media_neg, np.less)

    x_new = []
    y_new = []

    for i in range(len(x[0])):
        if (wy_media[t,x[0][i]]>0.5):
            x_new.append(x[0][i])

    for i in range(len(y[0])):
        if (wy_media[t,y[0][i]]<-0.5):
            y_new.append(y[0][i])


    plt.clf()
    plt.title(r'$\gamma$ = {:}'.format(int(gamma_array[g])), fontsize='20')
    plt.plot(z[x_new], np.transpose(wy_media[t,x_new]), 'ro')
    plt.plot(z[y_new], np.transpose(wy_media[t,y_new]), 'bo')
    # plt.plot(z[x], np.transpose(wy_media[t,x]), 'ro')
    # plt.plot(z[y], np.transpose(wy_media[t,y]), 'bo')
    plt.plot(z, wy_media[t, :])
    plt.xlim([0,2*np.pi])
    # plt.ylim([-0.03,0.03])
    plt.xlabel('z', fontsize='30')
    plt.ylabel(r'<$\omega_y$>', fontsize='30')
    plt.grid(True)
    plt.xticks(fontsize='20')
    plt.yticks(fontsize='20')
    plt.pause(0.01)
    plt.show()

    cant_max = len(x_new)
    cant_min = len(y_new)
    long_onda_max[g] = 8*np.pi / cant_max
    long_onda_min[g] = 8*np.pi / cant_min

# en primer lugar grafico solo para la longitud de onda maxima
plt.clf()
plt.title('Longitud de Onda del Modo mas Inestable', fontsize='20')
plt.plot(gamma_array, long_onda_min, 'x')
# plt.ylim([-0.03,0.03])
plt.xlabel(r'$\gamma$', fontsize='30')
plt.ylabel(r'$\lambda$', fontsize='30')
plt.grid(True)
plt.xticks(fontsize='20')
plt.yticks(fontsize='20')
plt.pause(0.01)
plt.show()
