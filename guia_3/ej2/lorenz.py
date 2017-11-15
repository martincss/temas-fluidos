# -*- coding: utf-8 -*-
"""
Martín Carusso, Florencia Lazzari
Temas avanzados de dinámica de fluidos, 2c 2017
Guía 3

Programa para la solución numérica del sistema Lorenz '63

"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import ode
# Ejemplo de uso de la clase 'ode'
#://stackoverflow.com/questions/12926393/using-adaptive-step-sizes-with-scipy-integrate-ode


#%% Definimos funciones para la integración

def F(t, X, sigma, r, b):
    """
    Dado el array (X,Y,Z) y los parámetros, devuelve 
    el array con las derivadas de (X,Y,Z) según las 
    ecuaciones de Lorenz 63
    """
    #X = np.array(X)
    DX = np.zeros(len(X))

    DX[0] = sigma*(X[1] - X[0])
    DX[1] = -X[0]*X[2] + r*X[0] - X[1]
    DX[2] = X[0]*X[1] - b*X[2]
    
    return DX


# Esto es una instance de la class ode, con el 
# parametro dopri5 en el method set_integrator
solver = ode(F).set_integrator('dopri5')

"""
El solver no requiere que se especifique un paso, se puede usar un callable
para que devuelva los valores a cada paso intermedio del integrate, como 
hicimos abajo
"""
# acá van a venir las soluciones a distintos
# tiempos
#list_sol = []

# necesitamos un callable para que reciba las
# soluciones del solver.integrate
#def call(t, x):
#    list_sol.append([t, *x])

#solver.set_solout(call)

# Se pasan los parámetros, la condición inicial y se integra el sistema

# Pasamos la lista de soluciones a un array
#list_sol = np.array(list_sol)
"""
En lugar de usar el callable para recibir la solución a cada paso intermedio, 
para poder controlar 'manualmente' la densidad de puntos, podemos pedir un
return de solver.integrate(t), para cada t de un array. Claramente es más 
lento, pero a nuestros efectos no es tan relevante como poder elegir la 
densidad de puntos.
"""

# meter el sigma, b y r
sigma = 10 
r = 30
b = 8/3
params = [sigma, r, b]

t0 = 0 # inicial
t1 = 50 # final, parece que con este solver no hay que poner el paso
W_0 = np.array([0, 0.5, 0.5])



def Lorenz(solver, params, W_0, tiempos):
    """
    Dado el solver para F, la lista de parámetros (sigma, r, b), la condición 
    inicial W_0, y un array de tiempos; integra el sistema Lorenz '63 para 
    cada tiempo del array.
    """
    # usamos methods en la clase ode para poner 
    # los valores iniciales y parámetros en la
    # función derivada
    solver.set_initial_value(W_0, tiempos[0])
    solver.set_f_params(*params)
    
    sol = np.zeros((len(tiempos), 3))
    for i in range(len(tiempos)):
        sol[i, :] = solver.integrate(tiempos[i])
    return sol
   
tiempos = np.linspace(t0, t1, 10000)
sol = Lorenz(solver, params, W_0, tiempos)

ax3 = Axes3D(plt.figure())
ax3.plot(sol[:,0], sol[:,1], sol[:,2])    
 
