import numpy as np
from cinpy.types import CDataTensor
from interp.engines import IntCub1DSet
import pyrho as rho
import os
import json
import logging
logging.basicConfig(level=logging.INFO)

exdir = os.path.abspath(os.path.join(__file__, os.pardir))
probefile = os.path.join(exdir,'probe.json')
datafile = os.path.join(exdir, 'data.bin')

# load in-vivo data
with open(probefile, 'r') as fp:
    probe = json.load(fp)

nang = int(11)
nele = int(probe['noElements'])
dele = probe['pitch']

alphas = (2*np.pi/180)*(np.arange(nang) - (nang-1)/2)

fc = 5208333.333 # Hz
fs = 4*fc
tstart = 5/fc
c = 1540
dt_lens = (probe['Rfocus'] - np.sqrt(probe["Rfocus"]**2 - probe["height"]**2/4))/c
dt_imp = 2/fc

# define the position of the elements
xele = dele*(np.arange(nele)-(nele-1)/2)
yele = np.zeros(nele)
zele = np.zeros(nele)
eles = np.array([xele, yele, zele])
eles = np.expand_dims(eles, axis=(1,3))

# calculate origins
origins = np.array([-dele*np.sign(alphas)*(nele-1/2), np.zeros(nang), np.zeros(nang)])
origins = np.expand_dims(origins, axis=(2,3))
origins = np.repeat(origins, axis=2, repeats=nele)

# define the reconstruction points in the field
xfield = np.linspace(-dele*(nele-1)/2, dele*(nele-1)/2, nele)
yfield = 0
zfield = np.arange(c*tstart, 40E-3, c/(2*fc))
Xfield, Yfield, Zfield = np.meshgrid(xfield, yfield, zfield, indexing='ij')
field = np.array([Xfield.flatten(), Yfield.flatten(), Zfield.flatten()])
field = np.expand_dims(field, axis=(1, 2))

# calculate lens offset
lens = np.expand_dims(np.array([0, 0, c*dt_lens]), axis=(1, 2, 3))

# calculate the normal vectors
norms = np.array([np.sin(alphas), np.zeros(nang), np.cos(alphas)])
norms = np.expand_dims(norms, axis=(2,3))

# calculate the TX tabs
tx = np.sum((field - origins + lens) * norms, axis=0)

# calculate the RX tabs
rx = np.linalg.norm((field - eles + lens), axis=0)

tau = rx+tx
del rx, tx, field, eles, lens, norms

print(tau.shape)

# load data
data_np = np.fromfile(datafile, dtype=np.int16).reshape((nang, nele, -1))
data = CDataTensor.fromnumpy(data_np)
print(data)
fint = IntCub1DSet(data, dn=1/fs, n0=tstart)

Nt = int(2*16+1)
points = np.array([0, *((np.arange(Nt)-(Nt-1)/2)/(8*fc))])

# beamformed = fint(tau, points)

# beamformed.copy2np().tofile("alldata.bin")
# print(beamformed.shape)