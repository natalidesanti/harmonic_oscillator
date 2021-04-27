'''
This program plots the comparison: classical probability X quantum density of probability in order to prove the Correspondence Principle
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os, shutil

files = ['fort.{}'.format(10+i) for i in range(8, 10)]
for f in files:
    shutil.copy(f, 'analytical/')
    os.remove(f)

files = ['fort.{}'.format(200+i) for i in range(31)]
for f in files:
    shutil.copy(f, 'analytical/')
    os.remove(f)

files = ['fort.{}'.format(600+i) for i in range(31)]
for f in files:
    shutil.copy(f, 'analytical/')
    os.remove(f)

data = {}
x = {}
y = {}
new_x = {}
new_y = {}
for i in range(30):
    data[i] = np.loadtxt(f'analytical/fort.{200 + i}') #Quantum
    x[i] = data[i][:,0]
    y[i] = data[i][:,1]
    new_x[i], new_y[i] = zip(*sorted(zip(x[i], y[i])))
ndata = {}
nx = {}
ny = {}
nnew_x = {}
nnew_y = {}
for i in range(30):
    ndata[i] = np.loadtxt(f'analytical/fort.{600 + i}') #Classical
    nx[i] = ndata[i][:,0]
    ny[i] = ndata[i][:,1]
    nnew_x[i], nnew_y[i] = zip(*sorted(zip(nx[i], ny[i])))

fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(nnew_x[0], nnew_y[0], label = f'Classical: $n = {0}$')
axs[0, 0].plot(new_x[0], new_y[0], label = f'Quantum: $n = {0}$')
axs[0, 0].set_xlim(-1.1, 1.1)
axs[0, 0].set_ylim(0., 3.0)
axs[0, 0].legend()
#
axs[0, 1].plot(nnew_x[1], nnew_y[1], label = f'Classical: $n = {1}$')
axs[0, 1].plot(new_x[1], new_y[1], label = f'Quantum: $n = {1}$')
axs[0, 1].set_xlim(-1.9, 1.9)
axs[0, 1].set_ylim(0., 3.0)
axs[0, 1].legend()
#
axs[1, 0].plot(nnew_x[2], nnew_y[2], label = f'Classical: $n = {2}$')
axs[1, 0].plot(new_x[2], new_y[2], label = f'Quantum: $n = {2}$')
axs[1, 0].set_xlim(-2.5, 2.5)
axs[1, 0].set_ylim(0., 3.0)
axs[1, 0].legend()
#
axs[1, 1].plot(nnew_x[3], nnew_y[3], label = f'Classical: $n = {3}$')
axs[1, 1].plot(new_x[3], new_y[3], label = f'Quantum: $n = {3}$')
axs[1, 1].set_xlim(-3., 3.)
axs[1, 1].set_ylim(0., 3.0)
axs[1, 1].legend()
#
for ax in axs.flat:
    ax.set(xlabel= r'$x$', ylabel= r'$|\psi_{n} (x)|^2, \rho_{n} (x)$')
#
for ax in axs.flat:
    ax.label_outer()
plt.savefig('analytical/correspondence_principle-0-3.png')
    
plt.figure(dpi = 100)
plt.plot(nnew_x[29], nnew_y[29], label = f'Classical: $n = {30}$')
plt.plot(new_x[29], new_y[29], label = f'Quantum: $n = {30}$')
plt.xlim(-8.5, 8.5)
plt.ylim(0., 1.0)
plt.legend()
plt.ylabel(r'$|\psi_{30} (x)|^2, \rho_{30} (x)$')
plt.xlabel(r'$x$')
plt.savefig('analytical/correspondence_principle-30.png')
