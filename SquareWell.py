# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 15:37:12 2013

@author: Nick Walker
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as an
import sys

i = 1j
l = 1.0

A = 0.0
B = 1.0

dx = 1e-3
x = np.linspace(0, l, 1000)

e = raw_input("range of energy eigenstates (space delimited, upper bound non-inclusive): ")

e = e.split()

w = (int(e[0]), int(e[1]))

phi = 0 * x

for n in range(w[0], w[1]):
    phi += np.sin(n * np.pi * x)

fig = plt.figure()
ax = plt.axes(xlim = (0, l), ylim = (-6, 6))

lines = [plt.plot([], [], lw = 2, color = 'blue'), plt.plot([], [], lw = 2, color = 'red'), plt.plot([], [], lw = 2, color = 'black')]

def init():
    for line, in lines:
        line.set_data([], [])
    return lines

def animate(k):
    t = 0.001 * k
    u = np.zeros((1, len(x)), dtype = float)
    v = np.zeros((1, len(x)), dtype = float)
    p = np.zeros((1, len(x)), dtype = float)
    s = [u, v, p]
    for n in range(w[0], w[1]):
        b = n * np.pi / l
        Cn = A * np.cos(b * x)
        Sn = B * np.sin(b * x)
        T = np.exp(-i * b ** 2 * t)
        An = 2.0 / l * sum(phi * Cn) * dx
        Bn = 2.0 / l * sum(phi * Sn) * dx
        u[0] = u[0] + np.real(T) * (An * Cn + Bn * Sn)
        v[0] = v[0] + np.imag(T) * (An * Cn + Bn * Sn)
    p[0] = u[0] ** 2 + v[0] ** 2
    n = 0
    for line, in lines:
        line.set_data(x, s[n][0] / (sum(p[0]) * dx))
        n += 1
    return lines

anim = an.FuncAnimation(fig, animate, init_func = init, frames = 10000, interval = 1)
plt.show()
#anim.save('square_wave_' + e[0] + '_' + e[1] + '.avi', writer = 'ffmpeg', fps = 60)
