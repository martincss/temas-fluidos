# -*- coding: utf-8 -*-
"""
Martín Carusso, Florencia Lazzari
Temas avanzados de dinámica de fluidos, 2c 2017
Guía 4

Programa para el análisis de los datos extraídos del código GHOST usando solver
HD (ejercicio 4)
"""
from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

f = open('balance.txt', 'r')

# delta_t_sampleo = 0.7
# tiempo = np.arrange(0, 20+delta_t_sampleo, delta_t_sampleo)
tiempo = np.zeros(154)
v2 = np.zeros(154)
w2 = np.zeros(154)
i = 0
for line in f:
    tiempo[i] = float(line[0:13])
    v2[i] = float(line[14:38])
    w2[i] = float(line[39:63])
    i += 1

plt.clf()
plt.title('Enstrofia en funcion del tiempo', fontsize='30')
plt.xlabel('Tiempo[s]', fontsize='30')
plt.ylabel(r'$\Omega$', fontsize='30')
plt.grid(True)
plt.plot(tiempo, w2)
plt.xlim([0,np.max(tiempo)])
plt.xticks(fontsize='20')
plt.yticks(fontsize='20')
plt.show()

# plt.plot(tiempo, v2)
# plt.show()

# calculo el maximo de enstrofia para obtener la tasa de disipacion
w2_max = np.max(w2)
t = np.where(w2 == w2_max)[0][0]

nu = 3.5e-3
disipacion = 2 * nu * w2_max

kolmogorov = (disipacion/nu**3)**(1./4)

onda_max_simulacion = 128./3

print 'Numero de onda de Kolmogorov: ',kolmogorov
print 'Maximo numero de onda resuelto por la simulacion: ',onda_max_simulacion

# promedio sobre 5 tiempos alrededor de t=5 (tiempo para el cual tenemos
# el maximo de entrofia)

# time step
delta_t = 7.0e-3
# number of steps between spectrum output
sstep = 100
cant_k = 65
k = np.linspace(1, 65, cant_k)

delta_spectrum = delta_t * sstep
t_spectrum = int(tiempo[t]/delta_spectrum)

tiempos = np.arange(t_spectrum-2, t_spectrum+3)
E_tiempos = np.zeros([len(tiempos), len(k)])

for j in range(len(tiempos)):
    t_spec = tiempos[j]
    f = open('kspectrum.000{:d}.txt'.format(t_spec), 'r')
    i = 0
    for line in f:
        E_tiempos[j, i] = float(line[13:36])
        i += 1

E_prom = np.mean(E_tiempos, axis=0)

# grafico espectro de kolmogorov (k**(-5/3))

kolmogorov = k**(-5./3)

plt.clf()
plt.title('Espectro de Energia Promediado', fontsize='30')
plt.xlabel('k', fontsize='30')
plt.ylabel('E', fontsize='30')
plt.grid(True)
plt.loglog(k, E_prom, label= 'Energia(k)')
plt.loglog(k[1:13], kolmogorov[1:13], label= r'$k^{-5/3}$')
plt.legend()
plt.xticks(fontsize='20')
plt.yticks(fontsize='20')
plt.show()

# estudio la funcion de transferencia

T_tiempos = np.zeros([len(tiempos), len(k)])

for j in range(len(tiempos)):
    t_spec = tiempos[j]
    f = open('ktransfer.000{:d}.txt'.format(t_spec), 'r')
    i = 0
    for line in f:
        T_tiempos[j, i] = float(line[13:36])
        i += 1

T_prom = np.mean(T_tiempos, axis=0)

plt.clf()
plt.title('Funcion Transferencia', fontsize='30')
plt.xlabel('k', fontsize='30')
plt.ylabel('T', fontsize='30')
plt.grid(True)
plt.plot(k, T_prom)
plt.xticks(fontsize='20')
plt.yticks(fontsize='20')
plt.show()

# calculo el flujo de energia

flujo_E = np.zeros(len(k))
vector_K = np.zeros(len(k))
suma = 0
suma_K = 0
for i in range(len(k)):
    suma_K += k[i]
    suma += T_prom[i]
    flujo_E[i] = -suma
    vector_K[i] = suma_K

plt.clf()
plt.title('Flujo de Energia', fontsize='30')
plt.xlabel('K', fontsize='30')
plt.ylabel(r'$\pi$', fontsize='30')
plt.grid(True)
plt.semilogy(vector_K, flujo_E)
plt.xticks(fontsize='20')
plt.yticks(fontsize='20')
plt.show()

# calculo E(t)

tiempo2 = 7*tiempo**-2

plt.clf()
plt.title('Energia en funcion del tiempo', fontsize='30')
plt.xlabel('Tiempo[s]', fontsize='30')
plt.ylabel('E', fontsize='30')
plt.grid(True)
plt.loglog(tiempo, v2/2, label= 'Energia(t)')
plt.loglog(tiempo[32:],tiempo2[32:], label=r'$t^{-2}$')
plt.xlim([0,np.max(tiempo)])
plt.legend()
plt.xticks(fontsize='20')
plt.yticks(fontsize='20')
plt.show()
