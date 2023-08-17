"""Python wrapper for c-type coherence functions"""
import ctypes as ct
import cinpy.types as cnpt
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

__trig__.pwtxengine.argtypes = ct.c_int, ct.c_float, ct.POINTER(ct.POINTER(ct.c_float)), ct.POINTER(ct.POINTER(ct.c_float)), ct.POINTER(ct.POINTER(ct.POINTER(ct.c_float))),
__trig__.pwtxengine.restype = ct.POINTER(ct.c_float)

__trig__.genmask3D.argtypes = ct.c_int, ct.c_float, ct.POINTER(ct.POINTER(ct.c_float)), ct.POINTER(ct.POINTER(ct.c_float)), ct.POINTER(ct.POINTER(ct.c_float)), ct.POINTER(ct.POINTER(ct.POINTER(ct.c_float))),
__trig__.genmask3D.restype = ct.POINTER(ct.c_int)

def geteletaus(eles, points, c:float=1540):
    """Calculate the time from each point to each element
    
    Parameters:
    ----
    eles: E by 3 matrix where E is the number of elements
    points: P by 3 matrix where P is the number of points
    c: speed of sound, default is 1540 M/s

    Returns:
    ----
    Taus: a list of pointers to E arrays, each holing the delay from the respective element to all the points
    """

    if (np.ndim(eles) != 2): raise ValueError("eles must be 2D")
    if (eles.shape[1] != 3):raise ValueError("eles must be E by 3")
    if (np.ndim(points) != 2): raise ValueError("points must be 2D")
    if (points.shape[1] != 3):raise ValueError("points must be E by 3")

    # extract constants
    # convert eles and points to c arrays
    Celes, Me, Ne = cnpt.copy2c(eles)
    Cpoints, Mp, Np = cnpt.copy2c(points)

    TauEs = []
    for i in range(Me.value):
        ptaue = __trig__.rxengine(Mp, ct.c_float(c), ct.byref(Celes[i]), ct.byref(Cpoints))
        TauEs.append(ptaue)

    return TauEs

