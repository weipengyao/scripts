import numpy as np

x = np.linspace(0, 160, 1600)

n0 = 55./2.*1.1e21

pres = 20.   # start position of preplasma
prel = 50.   # length of preplasma
prew = 1.   # width of the target

rear_start = pres+prel+prew  # start position of rear plasma
rear_length = 25.   # length of rear plasma

def preplasma(x):
    """
    x      : float, position [array]
    y      : float, position [array], for now, transverse profile is neglected
    n_max  : float, maximum preplasma density [which is the target density]
    n_min  : float, minimum preplasma density
    start  : float, start position of preplasma
    length : float, length of preplasma
    plateau: float, width of the target
    """
    start = pres
    length = prel
    plateau = prew
    n_max = n0 #dene
    n_min = 0.01*1.1e21
    if x <= start:
        return 0.
    elif x <= (start+length):
        return n_max*np.exp(np.log(n_max/n_min)*((x-(start+length))/length))
    elif x <= (start+length+plateau):
        return 1.*n_max
    else:
        return 0.

def rear_plasma(x):
    """
    x      : float, position [array]
    y      : float, position [array], for now, transverse profile is neglected
    n_max  : float, maximum preplasma density [which is the target density]
    n_min  : float, minimum preplasma density
    start  : float, start position of preplasma
    length : float, length of preplasma
    plateau: float, width of the target
    """
    start = rear_start
    length = rear_length
    # plateau = 0.
    n_max = n0 #dene
    n_min = 0.01*1.1e21
    if x <= start:
        return 0.
    elif x <= (start+length):
        return n_max*np.exp(np.log(n_max/n_min)*((start-x)/length))
    # elif x <= (start+length+plateau):
    #     return 1.*n_max
    else:
        return 0.

density = np.array([preplasma(xi)+rear_plasma(xi) for xi in x])



n_max = n0
n_min = 0.01*1.1e21

def preplasma2(x):
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

def rear_plasma2(x):
    if x < rear_start:
        return 0.0
    elif x <= rear_start + rear_length:
        # ramp down: n_max -> n_min
        t = (x - rear_start) / rear_length  # 0..1
        return n_max * np.exp(np.log(n_min/n_max) * t)
    else:
        return 0.0

density2 = np.array([preplasma2(xi) + rear_plasma2(xi) for xi in x])

import matplotlib.pyplot as plt
plt.plot(x, density, '-')
plt.plot(x, density2, '--')
plt.yscale('log')
plt.show()
