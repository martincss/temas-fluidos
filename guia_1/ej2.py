# -*- coding: utf-8 -*-
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt

def CD(x, delta_x):
	"""
	Hace la diferencia central para cdc periodicas
	"""
	diff = (np.roll(x,-1) - np.roll(x,+1))/(2*delta_x)
	return diff
	
def CD2(x, delta_x):
	"""
	Hace la segunda derivada para cdc periodicas 
	""" 
	diff2 = (np.roll(x,+1) - 2*x + np.roll(x,-1))/(delta_x**2) 
	return diff2

def F(x, delta_x, nu):
	"""
	Genera la funcion que da la derivada para Runge-Kutta 
	"""
	f = nu*CD2(x,delta_x) - x*(CD(x,delta_x))
	return f

def RK2(x, delta_x, delta_t, nu):
	"""
	Hace la integracion temporal por Runge-Kutta y devuelve la funcion
	en el tiempo proximo
	"""
	x_star = x + delta_t/2 * F(x, delta_x, nu)
	x_adv = x + delta_t * F(x_star, delta_x, nu)
	return x_adv 

#%%

x_min = 0
x_max = 2*pi
delta_x = 0.05
X = np.arange(x_min, x_max, delta_x)

t_inicial = 0
t_final = 12.03
delta_t = 0.001
T = np.arange(t_inicial, t_final, delta_t)

nu_1 = 0.1
nu_2 = 0.01
nu_3 = 0.001

u_1 = np.zeros((len(T), len(X)))
u_2 = np.zeros((len(T), len(X)))
u_3 = np.zeros((len(T), len(X)))

u_1[0,:] = np.sin(X)
u_2[0,:] = np.sin(X)
u_3[0,:] = np.sin(X)

#%%

for i in range(len(T)-1):
	u_1[i+1,:] = RK2(u_1[i,:], delta_x, delta_t, nu_1)
	u_2[i+1,:] = RK2(u_2[i,:], delta_x, delta_t, nu_2)
	u_3[i+1,:] = RK2(u_3[i,:], delta_x, delta_t, nu_3)

#%%

plt.figure()
for i in range(0, len(T)-1, 10):
    plt.clf()
    plt.xlim([0,2*pi])    
    plt.ylim([-1,1])
    plt.plot(X, u_1[i,:], 'b', label = '$\\nu = 0.1$'+' t = {:.2f}'.format(T[i]))
    plt.plot(X, u_2[i,:], 'g', label = '$\\nu = 0.01$'+' t = {:.2f}'.format(T[i]))
    plt.plot(X, u_3[i,:], 'r', label = '$\\nu = 0.001$'+' t = {:.2f}'.format(T[i]))
    plt.grid(True)
    plt.legend()
    plt.xlabel('$x$', fontsize = 15)
    plt.ylabel('$u$', fontsize = 15)
    plt.title('Solucion a la ecuacion de Burgers')
    plt.pause(0.001)
    #plt.savefig('/home/florlazz/Desktop/temas-fluidos/guia_1/frames_ej2/frame_{:03d}.png'.format(i))
plt.show()




