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

path = './'
f = open('balance.txt', 'r')

v2 = np.zeros(818)
w2 = np.zeros(818)
i = 0
for line in f:
    # fila[i] = line
    v2[i] = float(line[14:38])
    w2[i] = float(line[39:63])
    i += 1

plt.plot(w2)
plt.show()
