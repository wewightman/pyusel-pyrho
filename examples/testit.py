from pyrho.trig import geteletaus
import numpy as np
import cinpy as cnp
import ctypes as ct
import matplotlib.pyplot as plt


c = 1540
f = 5.2E6
nele = int(128)
dele = c/f

xele = dele*(np.arange(nele) - (nele-1)/2)
yele = np.zeros(nele)
zele = np.zeros(nele)
eles = np.array([xele, yele, zele]).T

fovrange = 30E-3
xfield = np.linspace(-fovrange/2, fovrange/2, int(2*fovrange/dele))
yfield = 0
zfield = np.linspace(dele , dele+fovrange, int(4*fovrange/dele))
Xfield, Yfield, Zfield = np.meshgrid(xfield, yfield, zfield, indexing='ij')
field = np.array([Xfield.flatten(), Yfield.flatten(), Zfield.flatten()]).T
print(field.shape)
taus = geteletaus(eles, field)
print(len(taus))
print(taus[0])
tau0 = cnp.copy2py(taus[nele//2], ct.c_int(field.shape[0]))
print(tau0.shape)

plt.figure()
plt.imshow(1E6*tau0.reshape(Xfield.squeeze().shape, order='c'))
plt.colorbar()

bins = np.linspace(0, 30, 30)
plt.figure()
plt.hist(1E6*tau0, bins=bins)
plt.show()