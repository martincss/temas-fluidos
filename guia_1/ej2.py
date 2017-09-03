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
	f = nu*CD2(x,delta_x) - x*(CD(x,delta_x))
	return f

def RK2(x, delta_x, delta_t, nu):
	x_star = x + delta_t/2 * F(x, delta_x, nu)
	x_adv = x + delta_t * F(x_star, delta_x, nu)
	return x_adv 

x_min = 0
x_max = 2*pi
delta_x = 0.2
X = np.arange(x_min, x_max, delta_x)

t_inicial = 0
t_final = 5
delta_t = 0.15
T = np.arange(t_inicial, t_final, delta_t)

nu = 0.1

u = np.zeros((len(T), len(X)))

u[0,:] = np.sin(X)

for i in range(len(T)-1):
	u[i+1,:] = RK2(u[i,:], delta_x, delta_t, nu)

fig = plt.figure()
ax = plt.axes(ylim = (-1, 1))
for i in range(len(T)-1):
    ax.plot(X, u[i,:], 'b', label = 't = {:.2f}'.format(T[i]))
    plt.ylim(-1, 1)
    plt.grid(True)
    plt.legend()
    plt.xlabel('$x$', fontsize = 15)
    plt.ylabel('$u$', fontsize = 15)
    plt.title('Solución a la ecuación de Burgers')
    plt.axis('tight')
    plt.pause(0.1)
    ax.clear()
    
plt.grid(True)
plt.legend()
plt.show()




