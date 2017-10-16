# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 17:30:11 2017

@author: martin
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
iters = 20

# Spatial resolution
NX = 128
NY = 128
NZ = 64
shape = (NX,NY,NZ)
x = np.linspace(0, 2*np.pi, NX)
delta_x = 2*np.pi/NX

field = np.fromfile(path+'th.{:04d}.out'.format(i),dtype=np.float32).reshape(shape,order='F')


#%% Fitteo
x_fit = x[:20]
temp_fit = field[:20, NY//2, NZ//2]
temp_fit /= np.mean(temp_fit)
f = lambda x, A, lon, ph: A*np.sin(2*np.pi/lon*x+ph)

p_guess = [4, 0.7, -0.2]
popt, _ = curve_fit(f, x_fit, temp_fit, p0 = p_guess)
fit = f(x_fit, *popt)

plt.plot(x_fit, fit, 'r')
plt.plot(x_fit, temp_fit, 'b')
plt.show()




