import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 160, 1600)

n0 = 55./2.*1.1e21
pres, prel, prew = 20., 50., 1.
rear_start, rear_length = pres+prel+prew, 25.

n_max = n0
n_min = 0.01*1.1e21

def preplasma(x):
    if x < pres:
        return 0.0
    elif x <= pres + prel:
        # ramp up: n_min -> n_max
        t = (x - pres) / prel  # 0..1
        return n_min * np.exp(np.log(n_max/n_min) * t)
    elif x <= pres + prel + prew:
        # solid target plateau
        return n_max
    else:
        return 0.0

def rear_plasma(x):
    if x < rear_start:
        return 0.0
    elif x <= rear_start + rear_length:
        # ramp down: n_max -> n_min
        t = (x - rear_start) / rear_length  # 0..1
        return n_max * np.exp(np.log(n_min/n_max) * t)
    else:
        return 0.0

density = np.array([preplasma(xi) + rear_plasma(xi) for xi in x])

plt.plot(x, density)
plt.yscale('log')
plt.xlabel('x')
plt.ylabel('density')
plt.show()
