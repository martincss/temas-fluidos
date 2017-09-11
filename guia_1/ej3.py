# -*- coding: utf-8 -*-
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import scipy.fftpack as fft

def term_adv(uk,k):
	du_dx = np.real(fft.ifft(1j*k*uk))
	ux = np.real(fft.ifft(uk))
	prod = ux*du_dx
	return fft.fft(prod)

def F(u,k,nu):
	f = -term_adv(u,k) - nu*k**2*u
	return f
	
def Euler(u, k, delta_t, nu):
    u_adv = u + F(u,k,nu)*delta_t
    return u_adv

def RK2(u, k, delta_t, nu):
	u_star = u + delta_t/2 * F(u, k, nu)
	u_adv = u + delta_t * F(u_star, k, nu)
	return u_adv

def Orszag(k, uk, N):
    for i in range(len(k)):
        if k[i] > N/3:
            uk[i] = 0
    

#%%

N = 2**6

"""
el programa es inestable para ciertas grillas, anda con N = 2**6 y 
delta_t = 0.02, pero no para delta_t más altos 
"""
x_min = 0
x_max = 2*pi
delta_x = 2*pi/N
X = np.arange(x_min, x_max, delta_x)

t_inicial = 0
t_final = 5
delta_t = 0.02
T = np.arange(t_inicial, t_final, delta_t)

nu = 0.1

u = np.zeros((len(T), len(X)))
uk = np.zeros((len(T), len(X)), dtype = complex)

u[0,:] = np.sin(X)
uk[0,:] = fft.fft(u[0,:])
#uk[0,:] = fft.fftshift(uk[0,:])
uk0 = uk[0,:]
k = 2 * np.pi * fft.fftfreq(N, delta_x)
#k = fft.fftfreq(N, delta_x)
#k = fft.fftshift(k)

#%%

for t in range(len(T)-2): 
    uk[t+1,:] = RK2(uk[t,:], k, delta_t, nu)
    #uk[t+1,:] = Euler(uk[t,:], k, delta_t, nu)
    Orszag(k, uk[t+1,:], N)
    u[t+1,:] = np.real(fft.ifft(uk[t+1,:]))

#%%

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






