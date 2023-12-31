"""Python wrapper for c-type coherence functions"""
import ctypes as ct
import cinpy as cnp
import numpy as np
from glob import glob
import platform as _pltfm
import os

# determine installed path
dirpath = os.path.dirname(__file__)

# determine the OS and relative binary file
if _pltfm.uname()[0] == "Windows":
    res = glob(os.path.abspath(os.path.join(dirpath, "*.dll")))
    name = res[0]
elif _pltfm.uname()[0] == "Linux":
    res = glob(os.path.abspath(os.path.join(dirpath, "*.so")))
    name = res[0]
else:
    res = glob(os.path.abspath(os.path.join(dirpath, "*.dylib")))
    name = res[0]

# load the c library
__trig__ = ct.CDLL(name)

__trig__.rxengine.argtypes = ct.c_int, ct.c_float, ct.POINTER(ct.POINTER(ct.c_float)), ct.POINTER(ct.POINTER(ct.POINTER(ct.c_float))),
__trig__.rxengine.restype = ct.POINTER(ct.c_float)

__trig__.pwtxengine.argtypes = ct.c_int, ct.c_float, ct.c_float, ct.POINTER(ct.POINTER(ct.c_float)), ct.POINTER(ct.POINTER(ct.c_float)), ct.POINTER(ct.POINTER(ct.POINTER(ct.c_float))),
__trig__.pwtxengine.restype = ct.POINTER(ct.c_float)

__trig__.genmask3D.argtypes = ct.c_int, ct.c_float, ct.POINTER(ct.POINTER(ct.c_float)), ct.POINTER(ct.POINTER(ct.c_float)), ct.POINTER(ct.POINTER(ct.c_float)), ct.POINTER(ct.POINTER(ct.POINTER(ct.c_float))),
__trig__.genmask3D.restype = ct.POINTER(ct.c_int)

def geteletaus(eles, points, c:float=1540):
    """Calculate the time from each point to each element
    
    Parameters:
    ----
    eles: E by 3 matrix where E is the number of elements
    points: P by 3 matrix where P is the number of points
    c: speed of sound, default is 1540 m/s

    Returns:
    ----
    Taus: a list of pointers to E arrays, each holing the delay from the respective element to all the points
    """

    if (np.ndim(eles) != 2): raise ValueError("eles must be 2D")
    if (eles.shape[1] != 3):raise ValueError("eles must be E by 3")
    if (np.ndim(points) != 2): raise ValueError("points must be 2D")
    if (points.shape[1] != 3):raise ValueError("points must be P by 3")

    # extract constants
    # convert eles and points to c arrays
    Celes, Me, Ne = cnp.copy2c(eles)
    Cpoints, Mp, Np = cnp.copy2c(points)

    TauEs = []
    for i in range(Me.value):
        ptaue = __trig__.rxengine(Mp, ct.c_float(c), ct.byref(Celes[i]), ct.byref(Cpoints))
        TauEs.append(ptaue)

    return TauEs

def getpwtaus(eles, points, teles, thetas, c=1540):
    """Calculate the delays for each element transition to each point - assuming plane wave tx

    Parameters:
    ----
    eles: location of each element, E by 3 matrix
    points: each point being interrogated, P by 3 matrix
    teles: time delay of each element at transmission, length E vector
    thetas: steering angle of each plane wave, length E vector
    c: assumed speed of sound

    Return:
    ----
    TauTx: a list of times needed to reach a given point, length E list of length P vectors
    """

    # convert eles and points to c arrays
    Celes, Me, Ne = cnp.copy2c(eles)
    Cpoints, Mp, Np = cnp.copy2c(points)
    Cteles, Mt = cnp.copy2c(teles.flatten())
    Cc = ct.c_float(c)

    if (Me.value != Mt.value) or (Mt.value != thetas.size):
        raise ValueError("All variables correxponding to the elements must have the same number of units")
    
    TauTx = []
    for ie in range(Me.value):
        pnorm = (ct.c_float * 3)(ct.c_float(np.sin(thetas[ie])), ct.c_float(np.cos(thetas[ie])), ct.c_float(0))
        tau = __trig__.pwtxengine(Mp, Cc, teles[ie], ct.byref(Celes[ie]), ct.byref(pnorm), ct.byref(Cpoints))
        TauTx.append(tau)

    return TauTx

def c_pw_engine(xref, tref, norm, points, Np:int, c:float=1540, dtype=ct.c_float):
    """Calculate the time to travel to each point based on spatial reference, temporal reference, 
    and normal vector for Ntx transmisions to Np points

    Parameters:
    ----
    xrefs: a list of POINTER(dtype) to vector of length 3, represents spatial reference point of wave
    trefs: a list of 'dtype' elements, represents the temporal reference point of the wave
    norms: a list of POINTER(dtype) to vector of length 3, represents the normal vector of wave propagation
    points: a POINTER(POINTER(dtype)) to matrix with size cNp 3 representing point in the field
    Np: number of points in field
    c: speed of sound in m/s
    dtype: ctype datatype - only currently supported option is ctype.c_float

    Return:
    ----
    TauTx: length Ntx list of POINTER(ctype) to vector of length Np. aka List of delaytabs to every point for each acq
    """

    cNp = ct.c_int(Np)
    cc = ct.c_float(c)

    ptaus = __trig__.pwtxengine(cNp, cc, tref, ct.byref(xref), ct.byref(norm), ct.byref(points))

    return ptaus

def c_norm_engine(xref, points, Np:int, c:float=1540, dtype=ct.c_float):
    """Calculate the time to travel from each point to each spatial reference

    Parameters:
    ----
    xref: a POINTER(dtype) to vector of length 3, represents spatial reference point of wave
    points: a POINTER(POINTER(dtype)) to matrix with size cNp 3 representing point in the field
    Np: number of points in field
    c: speed of sound in m/s
    dtype: ctype datatype - only currently supported option is ctype.c_float

    Return:
    ----
    TauTx: a POINTER(ctype) to vector of length Np. aka delaytabs to every point for each acq
    """

    cNp = ct.c_int(Np)
    cc = ct.c_float(c)

    ptaus = __trig__.rxengine(cNp, cc, ct.byref(xref), ct.byref(points))

    return ptaus
 
