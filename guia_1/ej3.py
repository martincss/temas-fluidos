# -*- coding: utf-8 -*-
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import scipy.fftpack as fft

def term_adv(uk,k):
	du_dx = fft.ifft(1j*k*uk)
	ux = fft.ifft(uk)
	prod = ux*du_dx
	return fft.fft(prod)

def F(u,k,nu):
	f = -term_adv(u,k) - nu*k**2*u
	return f
	

def RK2(u, k, delta_t, nu):
	u_star = u + delta_t/2 * F(u, k, nu)
	u_adv = u + delta_t * F(u_star, k, nu)
	return u_adv

N = 2**12

x_min = 0
x_max = 2*pi
X = np.linspace(x_min, x_max, N)
delta_x = X[1] - X[0]

t_inicial = 0
t_final = 5
delta_t = 0.15
T = np.arange(t_inicial, t_final, delta_t)

nu = 0.1

u = np.zeros((len(T), len(X)))
uk = np.zeros((len(T), len(X)), dtype = complex)

u[0,:] = np.sin(X)
uk[0,:] = fft.fft(u[0,:])
uk[0,:] = fft.fftshift(uk[0,:])
k = 2 * np.pi * fft.fftfreq(N, delta_x)
k = fft.fftshift(k)

for t in range(len(T)-2): 
    uk[t+1,:] = RK2(uk[t,:], k, delta_t, nu)
    u[t+1,:] = fft.ifft(uk[t,:])

plt.figure()
for i in range(len(T)-1):
    plt.clf()
    plt.xlim([0,2*pi])    
    plt.ylim([-1,1])
    plt.plot(X, u[i,:], 'b', label = 't = {:.2f}'.format(T[i]))
    plt.grid(True)
    plt.legend()
    plt.xlabel('$x$', fontsize = 15)
    plt.ylabel('$u$', fontsize = 15)
    plt.title('Solucion a la ecuacion de Burgers')
    plt.pause(0.001)
    #plt.savefig('/home/florlazz/Desktop/temas-fluidos/guia_1/frames_ej2/frame_{:03d}.png'.format(i))
plt.show()






