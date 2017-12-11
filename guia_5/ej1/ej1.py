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

#%%

path = 'data\\balance.txt'
datos = np.loadtxt(path)

