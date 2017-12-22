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
from scipy.integrate import simps

#%% Importamos los datos de balance

#path = 'data\\'
path = 'data/'
file_balance = 'balance.txt'
tiempo, energia, enstrofia, epsilon = np.loadtxt(path + file_balance, unpack = True)
#energia = 0.5*energia
#enstrofia = 0.5*enstrofia

#%% Graficamos estos datos

plt.figure()
plt.subplot(3,1,1)
plt.plot(tiempo, energia, color = 'r', label = 'Energía $\\langle v^2 \\rangle$')
#plt.xlabel('Tiempo')
plt.xticks([0, 5, 10, 15, 20, 25], ()) # para que no me muestre el label en los
plt.ylabel('Energía')                  # de arriba
plt.grid(True)
plt.legend()

plt.subplot(3,1,2)
plt.plot(tiempo, enstrofia, color = 'b', label = 'Enstrofía $\\langle \\omega^2 \\rangle$')
#plt.xlabel('Tiempo')
plt.xticks([0, 5, 10, 15, 20, 25], ())
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


#%%============================================================================
# ESPECTROS (punto c)
# =============================================================================

file_para = 'kspecpara.{:04d}.txt'
# desde el tiempo t* (ver el índice) promediamos los datos del espectro
index_star = 18
index_max = 50
cant_tiempos = index_max - index_star

k_para = np.loadtxt(path + file_para.format(1), delimiter = '  ', usecols = (0,))
E_para = np.zeros(len(k_para))

# desde el t* hasta el final, sumo todos los E(k) y después divido por cantidad
# de tiempos, así obtengo el promedio
for i in range(index_star, index_max):
    E_para += np.loadtxt(path + file_para.format(i), delimiter = '  ', usecols = (1,))
E_para /= cant_tiempos


# ahora para el perp
file_perp = 'kspecperp.{:04d}.txt'

k_perp = np.loadtxt(path + file_perp.format(1), delimiter = '  ', usecols = (0,))
E_perp = np.zeros(len(k_perp))

for i in range(index_star, index_max):
    E_perp += np.loadtxt(path + file_perp.format(i), delimiter = '  ', usecols = (1,))
E_perp /= cant_tiempos


# y el isótropo
file_iso = 'kspectrum.{:04d}.txt'

k_iso = np.loadtxt(path + file_iso.format(1), delimiter = '  ', usecols = (0,))
E_iso = np.zeros(len(k_iso))

for i in range(index_star, index_max):
    E_iso += np.loadtxt(path + file_iso.format(i), delimiter = '  ', usecols = (1,))
E_iso /= cant_tiempos

#%% Graficamos todo

plt.figure()
plt.subplot(1,3,1)
plt.loglog(k_iso, E_iso, color = 'b')
plt.xlabel('$k$')
plt.ylabel('$E(k)$')
plt.title('Espectro de energía isótropo')
plt.grid(True)
plt.legend()

plt.subplot(1,3,2)
plt.loglog(k_para, E_para, color = 'g')
plt.xlabel('$k_{\\parallel}$')
plt.ylabel('$E(k_{\\parallel})$')
plt.title('Espectro en la dirección paralela')
plt.grid(True)
plt.legend()

plt.subplot(1,3,3)
plt.loglog(k_perp, E_perp, color = 'r')
plt.xlabel('$k_{\\perp}$')
plt.ylabel('$E(k_{\\perp})$')
plt.title('Espectro en la dirección perpendicular')
plt.grid(True)
plt.legend()

#plt.tight_layout()


# =============================================================================
#%% LONGITUDES INTEGRALES (punto d)
# =============================================================================

def long_int(k, spectr, E_tot):
    '''
    Dada la k, y el espectro, teniendo la energía total, calcula la longitud
    integral correspondiente. Esta dada por L = 2pi int E(K)/k dk / E_tot
    '''
    integ = simps(spectr/k, x=k)
    L = 2*np.pi*integ/E_tot
    
    return L
    
L_iso = np.zeros(len(tiempo))
L_para = np.zeros(len(tiempo))    
L_perp = np.zeros(len(tiempo))    
    
for i in range(2, len(tiempo)+1):
    E_iso = np.loadtxt(path + file_iso.format(i), delimiter = '  ', usecols = (1,))    
    E_tot = simps(E_iso, x=k_iso)
    
    E_para = np.loadtxt(path + file_para.format(i), delimiter = '  ', usecols = (1,))
    E_perp = np.loadtxt(path + file_perp.format(i), delimiter = '  ', usecols = (1,))
    
    L_iso[i-1] = long_int(k_iso, E_iso, E_tot)
    L_para[i-1] = long_int(k_para[1:], E_para[1:], E_tot)
    L_perp[i-1] = long_int(k_perp[1:], E_perp[1:], E_tot)

#%% Graficamos esto
    
plt.figure()
plt.plot(tiempo[1:], L_iso[1:], color = 'b', label = '$L_{iso}$')
plt.plot(tiempo[1:], L_para[1:], color = 'g', label = '$L_{\\parallel}$')
plt.plot(tiempo[1:], L_perp[1:], color = 'r', label = '$L_{\\perp}$')
plt.xlabel('Tiempo')
plt.ylabel('Longitud integral')
plt.title('Longitudes integrales isótropa, paralela y perpendicular')
plt.grid(True)
plt.legend()
    
    
#==============================================================================
#%% CORTE DE OMEGA_Z EN EL PLANO X-Z (punto f) 
#==============================================================================
    
Nx = 192
Ny = 192
Nz = 48
shape = (Nx,Ny,Nz)
tiempo_wz = 12
wz = np.fromfile(path+'wz.{:04d}.out'.format(tiempo_wz),dtype=np.float32).reshape(shape,order='F')

plt.clf()
plt.title('Vorticidad $\\omega_z$ tiempo t = {:.1f}'.format(tiempo_wz*25/17), fontsize='20')
plt.imshow(np.transpose(wz[:,Ny//2,:]), cmap = 'hot', vmin = -1.2, vmax = 2.0, origin = 'lower')
plt.xlabel('x', fontsize='20')
plt.ylabel('z', fontsize='20')
plt.xticks(fontsize='20')
plt.yticks(fontsize='20')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=20)
plt.pause(0.01)
plt.show()
    
    
    







