# -*- coding: utf-8 -*-
"""
Martín Carusso, Florencia Lazzari
Temas avanzados de dinámica de fluidos, 2c 2017
Guía 3

Programa para la solución numérica del sistema Lorenz 63

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

#%% Definimos funciones para la integración

def F(X, sigma, r, b):
    """
    Dado el array (X,Y,Z) y los parámetros, devuelve 
    el array con las derivadas de (X,Y,Z) según las 
    ecuaciones de Lorenz 63
    """
    DX = np.zeros(len(X))

    DX[0] = sigma*(X[1] - X[0])
    DX[1] = -X[0]*X[2] + r*X[0] - X[1]
    DX[2] = X[0]*X[1] - B*X[2]
    
    return DX

# Esto es una instance de la class ode, con el 
# parametro dopri5 en el method set_integrator
solver = ode(F).set_integrator('dopri5')

# acá van a venir las soluciones a distintos
# tiempos
list_sol = []

# necesitamos un callable para que reciba las
# soluciones del solver.integrate
def call(t, x):
    list_sol.append([t, *x])

solver.set_solout(call)
solver.set_initial_value()

# meter el sigma, b y r
params = [0, 0, 0]

#://stackoverflow.com/questions/12926393/using-adaptive-step-sizes-with-scipy-integrate-ode


def Lorenz(solver)
    pass




