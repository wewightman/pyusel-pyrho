from pyrho.trig import geteletaus
import numpy as np
c = 1540
f = 5.2E6
nele = int(128)
dele = c/f

xele = dele*(np.arange(nele) - (nele-1)/2)
yele = np.zeros(nele)
zele = np.zeros(nele)
eles = np.array([xele, yele, zele]).T

range = 30E-3
xfield = np.linspace(-range/2, range/2, int(2*range/dele))
yfield = 0
zfield = np.linspace(dele , dele+range, int(4*range/dele))

