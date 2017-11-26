# -*- coding: utf-8 -*-
"""
Martín Carusso, Florencia Lazzari
Temas avanzados de dinámica de fluidos, 2c 2017
Guía 3

Programa para graficar la condicion inicial del ejercicio 1
"""

import numpy as np
from numpy import pi
import matplotlib.pyplot as plt

def u(z, u_0, gamma):
    U = u_0 * ( np.tanh(gamma * (z-np.pi/2)) + np.tanh(gamma * (-z + (3*np.pi)/2)) - 1)
    return U

u_0 = 1
gamma = 10
z = np.arange(0,2*np.pi, 0.01)
U_sol = u(z, u_0, gamma)

plt.clf()
plt.title('Condicion Inicial', fontsize='15')
plt.plot(z, U_sol)
plt.ylim([-1.2, 1.2])
plt.xlim([0 , 2*np.pi])
plt.xlabel(r'$z$', fontsize='20')
plt.ylabel(r'$u$', fontsize='20')
plt.legend()
plt.grid()
plt.xticks(fontsize='15')
plt.yticks(fontsize='15')
plt.show()

def du_dz2(z, u_0, gamma):
    U = (u_0 * gamma * ( (1./np.cosh(gamma * (z-np.pi/2)))**2 - (1./np.cosh(gamma * (-z + (3*np.pi)/2)))**2))**2
    return U

DU_DZ2_sol = du_dz2(z, u_0, gamma)

plt.clf()
plt.title(r'$(du/dz)^2$', fontsize='15')
plt.plot(z, DU_DZ2_sol)
plt.xlim([0 , 2*np.pi])
plt.xlabel(r'$z$', fontsize='20')
plt.ylabel(r'$u$', fontsize='20')
plt.legend()
plt.grid()
plt.xticks(fontsize='15')
plt.yticks(fontsize='15')
plt.show()
