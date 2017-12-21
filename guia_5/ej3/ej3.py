# -*- coding: utf-8 -*-
"""
Martín Carusso, Florencia Lazzari
Temas avanzados de dinámica de fluidos, 2c 2017
Guía 5

(ejercicio 3)
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 17})
from scipy.interpolate import interp1d

###############################################################
# a) graficamos energia total y helicidad cruzada en funcion
# del tiempo

path = 'data/'
file_balance = 'balance.txt'
tiempo, energia, omega2, j2 = np.loadtxt(path + file_balance, unpack = True)

file_cross = 'cross.txt'
tiempo_2, helicidad, a2 = np.loadtxt(path + file_cross, unpack = True)

plt.figure()
plt.loglog(tiempo, energia, color = 'b', label= u'$E(t)$')
plt.loglog(tiempo_2, helicidad, color = 'g', label= u'$K(t)$')
plt.xlabel('Tiempo', fontsize='25')
plt.title(u'Energía y Helicidad en función del Tiempo')
plt.grid(True)
plt.legend()
plt.show()

def central(x, y):
    '''
    Calcula mediante diferencia finita central, la derivada
    dy/dx. Para los valores de los bordes, repite los más cercanos.
    '''
    dydx = np.zeros(len(x))
    for i in range(1,len(x)-1):
        dydx[i] = (y[i+1] - y[i-1])/(x[i+1] - x[i-1])
    dydx[0] = dydx[1] # para mantener continuidad aprox
    dydx[-1] = dydx[-2] # idem

    return dydx

dE_dt = central(tiempo, energia)
dK_dt = central(tiempo_2, helicidad)

mod_dE_dt = abs(dE_dt)
mod_dK_dt = abs(dK_dt)

plt.figure()
plt.loglog(tiempo, mod_dE_dt, color = 'b', label= u'$\\frac{dE}{dt}$')
plt.loglog(tiempo_2, mod_dK_dt, color = 'g', label= u'$\\frac{dK}{dt}$')
plt.xlabel('Tiempo', fontsize='25')
plt.title(u'Tasa de Decaimiento')
plt.grid(True)
plt.legend(fontsize='20')
plt.show()



###############################################################
# b) graficamos energia cinetica y magnetica en funcion
# del tiempo

path = 'data/'
file_energy = 'energy.txt'
tiempo, v2, b2 = np.loadtxt(path + file_energy, unpack = True)

plt.figure()
plt.plot(tiempo, v2, color = 'b', label= u'$\\langle v^2 \\rangle$')
plt.plot(tiempo, b2, color = 'g', label= u'$\\langle b^2 \\rangle$')
plt.xlabel('Tiempo', fontsize='25')
plt.ylabel(u'Energía', fontsize='25')
plt.title(u'Energía en Función del Tiempo')
plt.grid(True)
plt.legend()
plt.show()

###############################################################
# c) graficamos <v.b> / (<v^2> <b^2>)^0.5 en funcion
# del tiempo


cociente = helicidad/(v2*b2)**0.5

plt.figure()
plt.plot(tiempo, cociente, color = 'b', label= u'$\\frac{\\langle \\vec{v} \\cdot \\vec{b}\\rangle}{\\left(\\langle {v}^2 \\rangle\\langle {b}^2\\rangle\\right)^\\frac{1}{2}}$')
plt.xlabel('Tiempo', fontsize='25')
plt.title(u'Valor medio de la Proyección entre $\\vec{v}$ y $\\vec{b}$')
plt.grid(True)
plt.legend()
plt.show()
