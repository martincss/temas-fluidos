# -*- coding: utf-8 -*-
"""
Martín Carusso, Florencia Lazzari
Temas avanzados de dinámica de fluidos, 2c 2017
Guía 5

Programa para graficar la energía, enstrofía y tasa de inyección de energía
(ejercicio 1b)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 15})

#%% Importamos los datos de balance

path = 'data\\'
file = 'balance.txt'
tiempo, energia, enstrofia, epsilon = np.loadtxt(path + file, unpack = True)
#energia = 0.5*energia
#enstrofia = 0.5*enstrofia

#%% Graficamos estos datos

plt.figure()
plt.subplot(3,1,1)
plt.plot(tiempo, energia, color = 'r', label = 'Energía $\\langle v^2 \\rangle$')
plt.xlabel('Tiempo')
plt.ylabel('Energía')
plt.grid(True)
plt.legend()

plt.subplot(3,1,2)
plt.plot(tiempo, enstrofia, color = 'b', label = 'Enstrofía $\\langle \\omega^2 \\rangle$')
plt.xlabel('Tiempo')
plt.ylabel('Enstrofía')
plt.grid(True)
plt.legend()

plt.subplot(3,1,3)
plt.plot(tiempo, epsilon, color = 'g', label = 'Inyección $\\epsilon$')
plt.xlabel('Tiempo')
plt.ylabel('Tasa de inyección \n de energía')
plt.grid(True)
plt.legend()

#plt.tight_layout()
#%%
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

dEdt = central(tiempo, energia)
nu = 2e-3
RHS = epsilon - 2*nu*enstrofia
plt.figure()
plt.plot(tiempo, dEdt, color = 'g', ls = '--', label = '$\\frac{dE}{dt}$')
plt.plot(tiempo, RHS, color = 'r', ls = '-.', label = '$\\epsilon - 2\\nu Z$')
plt.legend()
plt.grid(True)








