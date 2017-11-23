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

#%% COMUN A TODOS LOS ITEMS

t0 = 0 # inicial
t1 = 50 # final
W_0 = np.array([0, 0.5, 0.5])

tiempos = np.linspace(t0, t1, 10000)


#%%~~~~~~~~~~~~~~~~~~ ITEM a) ~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# meter el sigma, b y r
sigma = 10 
r_a = 2
b = 8/3
params_a = [sigma, r_a, b]

sol_a = Lorenz(solver, params_a, W_0, tiempos)

plt.figure()
plt.subplot(2,1,1)
plt.plot(tiempos, sol_a[:,1])
#plt.xlabel('Tiempo', fontsize = 15)
plt.ylabel('Coordenada $Y$', fontsize = 15)
plt.title('Integración de Lorenz con $r = {:d}$'.format(r_a), fontsize = 15)
plt.grid(True)

plt.subplot(2,1,2)
plt.plot(tiempos, sol_a[:,2])
plt.xlabel('Tiempo', fontsize = 15)
plt.ylabel('Coordenada $Z$', fontsize = 15)
#plt.title('Integración de Lorenz con $r = {:d}$'.format(r), fontsize = 15)
plt.grid(True)

plt.figure()
plt.plot(sol_a[:,0], sol_a[:,2])
plt.xlabel('Coordenada $Y$', fontsize = 15)
plt.ylabel('Coordenada $Z$', fontsize = 15)
plt.title('Integración de Lorenz con $r = {:d}$'.format(r_a), fontsize = 15)
plt.grid(True)

ax3 = Axes3D(plt.figure())
ax3.plot(sol_a[:,0], sol_a[:,1], sol_a[:,2])
plt.xlabel('Coordenada $X$', fontsize = 15)    
plt.ylabel('Coordenada $Y$', fontsize = 15)
ax3.set_zlabel('Coordenada $Z$', fontsize = 15)
ax3.set_title('Integración de Lorenz con $r = {:d}$'.format(r_a), fontsize = 15)

#%%~~~~~~~~~~~~~~~~~~ ITEM b_1) ~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# meter el sigma, b y r
sigma = 10 
r_b1 = 10
b = 8/3
params_b1 = [sigma, r_b1, b]

sol_b1 = Lorenz(solver, params_b1, W_0, tiempos)

plt.figure()
plt.subplot(2,1,1)
plt.plot(tiempos, sol_b1[:,1])
#plt.xlabel('Tiempo', fontsize = 15)
plt.ylabel('Coordenada $Y$', fontsize = 15)
plt.title('Integración de Lorenz con $r = {:d}$'.format(r_b1), fontsize = 15)
plt.grid(True)

plt.subplot(2,1,2)
plt.plot(tiempos, sol_b1[:,2])
plt.xlabel('Tiempo', fontsize = 15)
plt.ylabel('Coordenada $Z$', fontsize = 15)
#plt.title('Integración de Lorenz con $r = {:d}$'.format(r), fontsize = 15)
plt.grid(True)

plt.figure()
plt.plot(sol_b1[:,0], sol_b1[:,2])
plt.xlabel('Coordenada $Y$', fontsize = 15)
plt.ylabel('Coordenada $Z$', fontsize = 15)
plt.title('Integración de Lorenz con $r = {:d}$'.format(r_b1), fontsize = 15)
plt.grid(True)

ax3 = Axes3D(plt.figure())
ax3.plot(sol_b1[:,0], sol_b1[:,1], sol_b1[:,2])
plt.xlabel('Coordenada $X$', fontsize = 15)    
plt.ylabel('Coordenada $Y$', fontsize = 15)
ax3.set_zlabel('Coordenada $Z$', fontsize = 15)
ax3.set_title('Integración de Lorenz con $r = {:d}$'.format(r_b1), fontsize = 15)


#%%~~~~~~~~~~~~~~~~~~ ITEM b_2) ~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# meter el sigma, b y r
sigma = 10 
r_b2 = 24
b = 8/3
params_b2 = [sigma, r_b2, b]

sol_b2 = Lorenz(solver, params_b2, W_0, tiempos)

plt.figure()
plt.subplot(2,1,1)
plt.plot(tiempos, sol_b2[:,1])
#plt.xlabel('Tiempo', fontsize = 15)
plt.ylabel('Coordenada $Y$', fontsize = 15)
plt.title('Integración de Lorenz con $r = {:d}$'.format(r_b2), fontsize = 15)
plt.grid(True)

plt.subplot(2,1,2)
plt.plot(tiempos, sol_b2[:,2])
plt.xlabel('Tiempo', fontsize = 15)
plt.ylabel('Coordenada $Z$', fontsize = 15)
#plt.title('Integración de Lorenz con $r = {:d}$'.format(r), fontsize = 15)
plt.grid(True)

plt.figure()
plt.plot(sol_b2[:,0], sol_b2[:,2])
plt.xlabel('Coordenada $Y$', fontsize = 15)
plt.ylabel('Coordenada $Z$', fontsize = 15)
plt.title('Integración de Lorenz con $r = {:d}$'.format(r_b2), fontsize = 15)
plt.grid(True)

ax3 = Axes3D(plt.figure())
ax3.plot(sol_b2[:,0], sol_b2[:,1], sol_b2[:,2])
plt.xlabel('Coordenada $X$', fontsize = 15)    
plt.ylabel('Coordenada $Y$', fontsize = 15)
ax3.set_zlabel('Coordenada $Z$', fontsize = 15)
ax3.set_title('Integración de Lorenz con $r = {:d}$'.format(r_b2), fontsize = 15)


#%%~~~~~~~~~~~~~~~~~~ ITEM c) ~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# meter el sigma, b y r
sigma = 10 
r_b = 25
b = 8/3
params_c = [sigma, r_b, b]

sol_c = Lorenz(solver, params_c, W_0, tiempos)

plt.figure()
plt.subplot(2,1,1)
plt.plot(tiempos, sol_c[:,1])
#plt.xlabel('Tiempo', fontsize = 15)
plt.ylabel('Coordenada $Y$', fontsize = 15)
plt.title('Integración de Lorenz con $r = {:d}$'.format(r_b), fontsize = 15)
plt.grid(True)

plt.subplot(2,1,2)
plt.plot(tiempos, sol_c[:,2])
plt.xlabel('Tiempo', fontsize = 15)
plt.ylabel('Coordenada $Z$', fontsize = 15)
#plt.title('Integración de Lorenz con $r = {:d}$'.format(r), fontsize = 15)
plt.grid(True)

plt.figure()
plt.plot(sol_c[:,0], sol_c[:,2])
plt.xlabel('Coordenada $Y$', fontsize = 15)
plt.ylabel('Coordenada $Z$', fontsize = 15)
plt.title('Integración de Lorenz con $r = {:d}$'.format(r_b), fontsize = 15)
plt.grid(True)

ax3 = Axes3D(plt.figure())
ax3.plot(sol_c[:,0], sol_c[:,1], sol_c[:,2])
plt.xlabel('Coordenada $X$', fontsize = 15)    
plt.ylabel('Coordenada $Y$', fontsize = 15)
ax3.set_zlabel('Coordenada $Z$', fontsize = 15)
ax3.set_title('Integración de Lorenz con $r = {:d}$'.format(r_b), fontsize = 15)



#%%~~~~~~~~~~~~~~~~~~ ITEM d) ~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# meter el sigma, b y r
sigma = 10 
r_b = 30
b = 8/3
params_d = [sigma, r_b, b]

sol_d = Lorenz(solver, params_d, W_0, tiempos)

W_0d = np.array([0, 0.5, 0.50001])
sol_dd = Lorenz(solver, params_d, W_0d, tiempos)


plt.figure()
plt.plot(tiempos, sol_d[:,1], label = '$z_0 = 0.5$')
plt.plot(tiempos, sol_dd[:,1], label = '$z_0 = 0.50001$')
plt.xlabel('Tiempo', fontsize = 15)
plt.ylabel('Coordenada $Y$', fontsize = 15)
plt.title('Integración de Lorenz con $r = {:d}$'.format(r_b), fontsize = 15)
plt.grid(True)
plt.legend()

ax3 = Axes3D(plt.figure())
ax3.plot(sol_d[:,0], sol_d[:,1], sol_d[:,2], label = '$z_0 = 0.5$', color = 'r')
ax3.plot(sol_dd[:,0], sol_dd[:,1], sol_dd[:,2], label = '$z_0 = 0.50001$', color = '#00B7EB')
plt.xlabel('Coordenada $X$', fontsize = 15)    
plt.ylabel('Coordenada $Y$', fontsize = 15)
ax3.set_zlabel('Coordenada $Z$', fontsize = 15)
ax3.set_title('Integración de Lorenz con $r = {:d}$'.format(r_b), fontsize = 15)
plt.legend()

    


