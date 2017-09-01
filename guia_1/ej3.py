# -*- coding: utf-8 -*-
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from scipy.fftpack import fft

def term_adv(uk,k):
	du_dx = fft.ffti(1i*k*uk)
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
delta_x = x[1] - x[0]

t_inicial = 0
t_final = 5
delta_t = 0.15
T = np.arange(t_inicial, t_final, delta_t)

nu = 0.1

u = np.zeros((len(T), len(X)))
uk = np.zeros((len(T), len(X)))

u[0,:] = np.sin(X)
u_k[0,:] = fft.fft(u[0,:])
k = 2 * np.pi * fft.fftfreq(N, delta_x)

for i in range(len(T)-2):
	for k in 
		uk[i+1,:] = RK2(uk[i,:], , delta_t, nu)
 






