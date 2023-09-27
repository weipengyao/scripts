#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 21:43:00 2023

@author: yao
"""

import numpy as np
import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt

plt.rc('font', family='sans-serif', serif='Arial')
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('axes', labelsize=12)

width  = 3.487
height = width / 1.618


xx = np.linspace(0,425,426) * 4.9 / 1e3
yy = np.loadtxt('/Users/yao/Desktop/test.txt') * 4.21e21 / 1e22


fig, ax = plt.subplots()
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)

ax.plot(xx, yy, color='k', ls='solid',lw=2.0)
ax.set_yticks([0, 0.5, 1.0, 1.5, 2.0])

ax.set_xlabel('x (mm)',fontsize=12)
ax.set_ylabel('$\int n_e dy$ ($10^{22} m^{-2}$)',fontsize=12)
ax.set_xlim([0,2])
ax.set_ylim([0,2.0])

fig.set_size_inches(width, height)
# plt.show()
fig.savefig('/Users/yao/Desktop/plot.pdf')