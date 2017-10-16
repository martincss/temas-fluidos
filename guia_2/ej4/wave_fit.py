# -*- coding: utf-8 -*-
"""
Martín Carusso, Florencia Lazzari
Temas avanzados de dinámica de fluidos, 2c 2017
Guía 2

Programa para la determinacion de la longitud de ondas de sotavento en función
de la frecuencia de Brunt-Väisälä a partir de soluciones a las ecuaciones
de Boussinesq (ejercicio 4b)
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.fftpack as fft

#%% Define el fitteo por un seno estimando la frecuencia según el pico de 
#   fourier

def fit_sin(tt, yy):
    """
    Ajusta seno a la señal y devuelve un diccionario con los parámetros del
    fit, y la función de ajuste.
    Crédito: https://stackoverflow.com/a/42322656
    """
    tt = np.array(tt)
    yy = np.array(yy)
    ff = fft.fftfreq(len(tt), (tt[1]-tt[0]))   # asume espaciamiento uniforme
    Fyy = abs(fft.fft(yy))
    # Para estimar la frecuencia, toma el pico del abs del espectro del fft
    guess_freq = abs(ff[np.argmax(Fyy[1:])+1])   # excluye el pico a frecuencia cero, relacionado con el offset
    guess_amp = (np.max(yy) - np.min(yy))*0.5
    guess_offset = np.mean(yy)
    guess = np.array([guess_amp, 2.*np.pi*guess_freq, 0., guess_offset])

    def sinfunc(t, A, w, p, c):  return A * np.sin(w*t + p) + c
    
    popt, pcov = curve_fit(sinfunc, tt, yy, p0=guess)
    A, w, p, c = popt
    f = w/(2.*np.pi)
    fitfunc = lambda t: A * np.sin(w*t + p) + c
    
    return {"amp": A, "omega": w, "phase": p, "offset": c, "freq": f, 
            "period": 1./f, "fitfunc": fitfunc, "maxcov": np.max(pcov), 
            "rawres": (guess,popt,pcov)}


#%% Preliminares para los archivos binarios

# Path a los directorios de cada BV
path = './N_{:d}/'
BVs = [1, 3, 5, 7, 10, 12, 15] # cuando estén todas cambiar por un range o enum
samples = 11                   # cantidad de samples para un N particular
longs = []                     # lista a contener las longitudes de onda

# Spatial resolution
NX = 128
NY = 128
NZ = 64
shape = (NX,NY,NZ)
x = np.linspace(0, 2*np.pi, NX)
delta_x = 2*np.pi/NX


#%% Fitteo para determinar lambda para cada N
"""
Para cada N, del directorio correspondiente toma el .out de la última iteración
en la cual se ve la onda ya formada. Toma el primer período de la onda y 
ajusta por un seno, en base a las estimaciones de los parámetros según la 
función fit_sin. Regista la longitud de onda a una lista.
"""

for N in BVs:
    field = np.fromfile(path.format(N) + 'th.{:04d}.out'.format(samples),
                        dtype=np.float32).reshape(shape,order='F')
    field = field[:, NY//2, NZ//2] # fijo Y,Z a la mitad del recinto
    # para el fit, me quedo sólo con un periodo de la señal (no me queda 
    # claro por qué funciona con *3 en lugar de *2)                                                         
    period = abs(np.argmax(field) - np.argmin(field))*3 
    # Si el periodo es más grande que el tamaño de la señal, tomo el más chico
    samp_fit = np.min([len(field), period])
    x_fit = x[:samp_fit]
    temp_fit = field[:samp_fit]
    # Normalizo la señal a la unidad, para mejorar el fit sobre la longitud 
    # de onda. Comentar la línea si se quiere la amplitud de verdad
    temp_fit /= abs(np.max(temp_fit) - np.min(temp_fit))*0.5
    res = fit_sin(x_fit, temp_fit)
    longs.append(res['period'])
    
    #Grafica el período de la señal con su fit, para verificar
    plt.plot(x_fit, temp_fit, 'r')
    plt.plot(x_fit, res['fitfunc'](x_fit), 'b')

#%% Plotteo de los resultados finales

N = np.array(BVs)
longs = np.array(longs)
plt.figure()
plt.plot(N, longs, 'g')
plt.title('Longitud de la onda estacionaria en función de N')
plt.xlabel('Frecuencia de Brunt-Väisälä')
plt.ylabel('Longitud de onda')
plt.grid(True)
















