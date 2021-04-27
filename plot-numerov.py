'''
This program plots the comparison: Numerov's X analytical probability densities
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os, shutil

shutil.copy('fort.666', 'numerov/')
os.remove('fort.666')
files = ['fort.{}'.format(700+i) for i in range(31)]
for f in files:
    shutil.copy(f, 'numerov/')
    os.remove(f)
files = ['fort.{}'.format(800+i) for i in range(31)]
for f in files:
    shutil.copy(f, 'numerov/')
    os.remove(f)

data = {}
x = {}
y = {}
new_x = {}
new_y = {}
for i in range(30):
    data[i] = np.loadtxt(f'analytical/fort.{200 + i}')
    x[i] = data[i][:,0]
    y[i] = data[i][:,1]
    new_x[i], new_y[i] = zip(*sorted(zip(x[i], y[i])))
ndata = {}
nx = {}
ny = {}
nnew_x = {}
nnew_y = {}
for i in range(30):
    ndata[i] = np.loadtxt(f'numerov/fort.{800 + i}')
    nx[i] = ndata[i][:,0]
    ny[i] = ndata[i][:,1]
    nnew_x[i], nnew_y[i] = zip(*sorted(zip(nx[i], ny[i])))

x = np.linspace(-2.5, 2.5, 100)
analytical = {}
numerov = {}
fig = plt.figure(constrained_layout=False, dpi=100)
gs = fig.add_gridspec(nrows=7, ncols=1)
f_ax1 = fig.add_subplot(gs[0:5,0])
for i in range(3):
    plt.plot(nnew_x[i], nnew_y[i], label = f'Numerov: $n = {i}$')
    plt.plot(new_x[i], new_y[i], '--', label = f'Analytical: $n = {i}$')
plt.xlim(-2.5, 2.5)
plt.ylim(0, 0.6)
plt.legend(loc = 'upper left')
plt.xlabel(r'$x$')
plt.ylabel(r'$|\psi_n (x)|^2$')
f_ax2 = fig.add_subplot(gs[5:7, 0])
for i in range(3):
    analytical[i] = interp1d(new_x[i], new_y[i])
    numerov[i] = interp1d(nnew_x[i], nnew_y[i])
    plt.plot(x, abs(1. - analytical[i](x)/numerov[i](x)), label = r'$\langle \Delta \rangle = {:.3f}$'.format(np.mean(abs(1. - analytical[i](x)/numerov[i](x)))))
plt.legend(loc = 'upper left')
plt.xlim(-2.5, 2.5)
plt.ylabel(r'Residue')
plt.xlabel(r'$x$')
plt.savefig('numerov/numerovXanalytical-0_2.png')

analytical30 = interp1d(new_x[29], new_y[29])
numerov30 = interp1d(nnew_x[29], nnew_y[29])
x = np.linspace(-8.5, 8.5, 100)

fig = plt.figure(constrained_layout=False, dpi=100)
gs = fig.add_gridspec(nrows=7, ncols=1)
f_ax1 = fig.add_subplot(gs[0:5,0])
plt.plot(nnew_x[29], nnew_y[29], label = f'Numerov: $n = {30}$')
plt.plot(new_x[29], new_y[29], '--', label = f'Analytical: $n = {30}$')
plt.xlim(-8.5, 8.5)
plt.legend()
plt.ylabel(r'$|\psi_n (x)|^2$')
f_ax2 = fig.add_subplot(gs[5:7, 0])
plt.plot(x, abs(1. - analytical30(x)/numerov30(x)), label = r'$\langle \Delta \rangle = {:.3f}$'.format(np.mean(abs(1. - analytical30(x)/numerov30(x)))))
plt.legend()
plt.ylabel(r'Residue')
plt.xlabel(r'$x$')
plt.savefig('numerov/numerovXanalytical-30.png')
