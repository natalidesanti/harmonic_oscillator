'''
This program plots the Hermite polynomials until order 5
'''

import numpy as np
import matplotlib.pyplot as plt
import os, shutil

files = ['fort.{}'.format(i+10) for i in range(6)]
for f in files:
    shutil.copy(f, 'hermite/')
    os.remove(f)

data = {}
for i in range(5):
    data[i] = np.loadtxt(f'hermite/fort.{i+10}')

plt.figure(dpi = 100)
plt.title('Hermite polynomials')
for i in range(5):
    plt.plot(data[i][:, 0], data[i][:, 1], label = f'$n = {i}$')
plt.legend()
plt.ylabel(r'$H_n (x)$')
plt.xlabel(r'$x$')
plt.xlim(-2.5, 2.5)
plt.ylim(-60, 60)
plt.savefig('hermite/plot_hermite.png')
