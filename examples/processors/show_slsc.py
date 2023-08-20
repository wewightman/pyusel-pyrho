import numpy as np
import matplotlib.pyplot as plt
from ctypes import c_float
from pyrho.processors.slsc import PWTX, FullRX, SLSCProc
import os
import json

exdir = os.path.abspath(os.path.join(__file__, os.pardir))
probefile = os.path.join(exdir,'probe.json')

# load in-vivo data
with open(probefile, 'r') as fp:
    probe = json.load(fp)

nang = int(11)
nele = int(probe['noElements'])
dele = probe['pitch']

fc = 5208333.333 # Hz
fs = 4*fc
tstart = 5/fc
c = 1540

# define the position of the elements
xele = dele*(np.arange(nele)-(nele-1)/2)
yele = np.zeros(nele)
zele = np.zeros(nele)
eles = np.array([xele, yele, zele]).T

# define the reconstruction points in the field
xfield = np.linspace(-dele*(nele-2/3)/2, dele*(nele-2/3)/2, 2*nele)
yfield = 0
zfield = np.arange(c*tstart, 40E-3, c/(2*fc))
Xfield, Yfield, Zfield = np.meshgrid(xfield, yfield, zfield, indexing='ij')
field = np.array([Xfield.flatten(), Yfield.flatten(), Zfield.flatten()]).T

# define the transmission/steering angles
alphas = 1.5*np.pi*(np.arange(nang)-(nang-1)/2)/180
xrefs = [eles[0,:] if alpha >= 0 else eles[-1,:] for alpha in alphas]
xrefs = np.array(xrefs)
trefs = np.zeros(nang)

tx = PWTX(alphas, xrefs, trefs, c=c, dtype=c_float)
rx = FullRX(eles, c=c, dtype=c_float)

loc = SLSCProc(c, field, tx, rx, lags=1)

loc(np.zeros((5,5)))