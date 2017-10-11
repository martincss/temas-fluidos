# -*- coding: utf-8 -*-
"""
Martín Carusso, Florencia Lazzari
Temas avanzados de dinámica de fluidos, 2c 2017
Guía 1

Programa para la solución numérica de la ecuación de Burgers por métodos
pseudo-espectrales (ejercicios 3,4)
"""
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import scipy.fftpack as fft

def term_adv(uk,k):
    """
    Calcula el término de advección en el espacio real y lo devuelve en el
    espacio de Fourier
    """
    du_dx = np.real(fft.ifft(1j*k*uk))
    ux = np.real(fft.ifft(uk))
    prod = ux*du_dx
    return fft.fft(prod)

def F(u,k,nu):
    """
    Genera la funcion que da la derivada para Runge-Kutta
    """
    f = -term_adv(u,k) - nu*k**2*u
    return f

def Euler(u, k, delta_t, nu):
    """
    Integración numérica por Euler (prueba para debugging de Runge-Kutta)
    """
    u_adv = u + F(u,k,nu)*delta_t
    return u_adv

def RK2(u, k, delta_t, nu):
    """
    Hace la integracion temporal por Runge-Kutta y devuelve la funcion
    en el tiempo proximo
    """
    u_star = u + delta_t/2 * F(u, k, nu)
    u_adv = u + delta_t * F(u_star, k, nu)
    return u_adv

def Orszag(k, uk, N):
    """
    Elimina las amplitudes de Fourier para los modos k > N/3
    (regla de los 2/3 de Orszag)
    """
    dealiasing = np.abs(k) < N/3
    uk *= dealiasing


#%% Inicializamos discretización espacial y temporal
####################################################

N = 2**8

x_min = 0
x_max = 2*pi
delta_x = 2*pi/N
X = np.arange(x_min, x_max, delta_x)

t_inicial = 0
t_final = 5
delta_t = 0.001
T = np.arange(t_inicial, t_final, delta_t)

nu = 0.1

u = np.zeros((len(T), len(X)))
uk = np.zeros((len(T), len(X)), dtype = complex)

# Condiciones iniciales
u[0,:] = np.sin(X)
k = 2 * np.pi * fft.fftfreq(N, delta_x)
uk[0,:] = fft.fft(u[0,:])
Orszag(k, uk[0,:], N)

#uk[0,:] = fft.fftshift(uk[0,:])
uk0 = uk[0,:]

#k = fft.fftfreq(N, delta_x)
#k = fft.fftshift(k)

#%% Integración temporal por Runge-Kutta de orden 2
####################################################

for t in range(len(T)-2):
    uk[t+1,:] = RK2(uk[t,:], k, delta_t, nu)
    #uk[t+1,:] = Euler(uk[t,:], k, delta_t, nu)
    Orszag(k, uk[t+1,:], N)
    u[t+1,:] = np.real(fft.ifft(uk[t+1,:]))

#%% Graficamos las soluciones
####################################################
"""
El último argumento del range cambia cada cuántas iteraciones temporales se
actualiza el gráfico. El plt.pause() permite graficar "en tiempo real". Para
guardar las frames, comentar plt.pause y descomentar la línea con plt.savefig
indicando un path acorde.
"""

#plt.figure()
#for i in range(0, len(T)-1, 10):
#    plt.clf()
#    plt.xlim([0,2*pi])
#    plt.ylim([-1,1])
#    plt.plot(X, u[i,:], 'b', label = 't = {:.2f}'.format(T[i]))
#    plt.grid(True)
#    plt.legend()
#    plt.xlabel('$x$', fontsize = 15)
#    plt.ylabel('$u$', fontsize = 15)
#    plt.title('Solucion a la ecuacion de Burgers')
#    plt.pause(0.001)
#    #plt.savefig('/home/florlazz/Desktop/temas-fluidos/guia_1/frames_ej2/frame_{:03d}.png'.format(i))
#plt.show()

"""
Dentro del loop se calcula la energia cinetica para cada valor de tiempo t, a
partir de la u(k) (transformada de Fourier de la u(x). El primer grafico es
energia cinetica en funcion de k. El segundo grafico tiene los mismos datos
graficados en modo loglog.
"""

k = fft.fftshift(k)

plt.figure()
for t in range(0, len(T)-1, 100):
   plt.clf()
   plt.ylim([0,9000])
   plt.xlim([-5,5])
   Ek = np.zeros(len(uk[t,:]))
   for i in range(len(uk[t,:])):
       Ek[i] = (np.absolute(uk[t,i]))**2 / 2
   Ek = fft.fftshift(Ek)
   plt.plot(k, Ek, 'b', label = 't = {:.2f}'.format(T[t]))
   plt.grid(True)
   plt.legend()
   plt.xlabel('$k$', fontsize = 15)
   plt.ylabel('$E$', fontsize = 15)
   plt.title('Energia en funcion de k para $\\nu$ = {:.2f}'.format(nu))
   plt.pause(0.001)
   #plt.savefig('/home/florlazz/Desktop/temas-fluidos/guia_1/frames_ej2/frame_{:03d}.png'.format(i))
plt.show()

# plt.figure()
# for t in range(0, len(T)-1, 100):
#    plt.clf()
#    plt.ylim([10**(-40),10**5])
#    Ek = np.zeros(len(uk[t,:]))
#    for i in range(len(uk[t,:])):
#        Ek[i] = (np.absolute(uk[t,i]))**2 / 2
#    Ek = fft.fftshift(Ek)
#    plt.loglog(k, Ek, 'b', label = 't = {:.2f}'.format(T[t]))
#    plt.grid(True)
#    plt.legend()
#    plt.xlabel('$k$', fontsize = 15)
#    plt.ylabel('$E$', fontsize = 15)
#    plt.title('Energia en funcion de k para $\\nu$ = {:.2f}'.format(nu))
#    plt.pause(0.001)
#    #plt.savefig('/home/florlazz/Desktop/temas-fluidos/guia_1/frames_ej2/frame_{:03d}.png'.format(i))
# plt.show()


"""
Se calcula la energia cinetica en funcion de u(x) segun la relacion E = u**2 / 2
y se grafica E en funcion de t.
"""

# E = np.zeros(len(T))
# for i in range(len(T)):
# 	for j in range(len(X)):
# 		E_x = u[i,j]**2 / 2
# 		E[i] += E_x

# plt.plot(E, T, 'b')
# plt.grid(True)
# plt.legend()
# plt.xlabel('$t$', fontsize = 15)
# plt.ylabel('$E$', fontsize = 15)
# plt.title('Energia en funcion de t para $\\nu$ = {:.2f}'.format(nu))
# plt.show()

"""
Genera 4 subplots donde se grafica la solucion de burguers, u en funcion de x
para distintos instantes temporales.
"""

# f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=True)
# plt.xlim([0,2*pi])
# plt.ylim([-1,1])
# ax1.plot(X, u[0,:], 'b', label = 't = {:.2f}'.format(T[0]))
# ax1.set_title('Solucion a la ecuacion de Burgers para $\\nu$ = 0.001')
# ax1.grid(True)
# ax1.legend()
# ax2.plot(X, u[40,:], 'b', label = 't = {:.2f}'.format(T[40]))
# ax2.grid(True)
# ax2.legend()
# ax2.set_ylabel('$u$', fontsize = 17)
# ax3.plot(X, u[100,:], 'b', label = 't = {:.2f}'.format(T[100]))
# ax3.grid(True)
# ax3.legend()
# ax4.plot(X, u[180,:], 'b', label = 't = {:.2f}'.format(T[180]))
# ax4.grid(True)
# ax4.legend()
# ax4.set_xlabel('$x$', fontsize = 17)
# # Fine-tune figure; make subplots close to each other and hide x ticks for
# # all but bottom plot.
# f.subplots_adjust(hspace=0)
# plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
# plt.legend()
# plt.show()
