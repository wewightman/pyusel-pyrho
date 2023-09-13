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
__rho__ = ct.CDLL(name)

__rho__.lagNRho.argtypes = ct.POINTER(ct.POINTER(ct.POINTER(ct.c_float))), ct.c_int, ct.c_int, ct.c_int
__rho__.lagNRho.restype = ct.c_float

# lagNRhoSet(int lag, int Np, int Ntx, int Nrx, int Nt, float Ts, float tstart, float dtw, float *** tautx, float *** taurx, float **** rf)
__rho__.lagNRhoSet.argtypes = (
    ct.c_int, 
    ct.c_int, 
    ct.c_int, 
    ct.c_int, 
    ct.c_int,
    ct.c_float,
    ct.c_float,
    ct.c_float,
    ct.POINTER(ct.POINTER(ct.POINTER(ct.c_float))),
    ct.POINTER(ct.POINTER(ct.POINTER(ct.c_float))),
    ct.POINTER(ct.POINTER(ct.POINTER(ct.POINTER(ct.c_float)))),)
__rho__.lagNRhoSet.restype = ct.POINTER(ct.c_float)

def lagNRho(arr, lag:int=1, axis:int=1):
    """calculate the Nth lag of 2D input matrix
    
    Parameters:
    ----
    arr: the array over which the Nth lag coherence is being calculated
    lag: the lag index
    axis: the axis corresponding to the elements

    Returns:
    ----
    rho: normalized correlation coefficient at this lag
    """

    # only allow 2D input arrays
    if np.ndim(arr) != 2:
        raise ValueError("Input array must be 2D")
    
    # Make the sample axis the fast axis
    if axis == 1:
        c_arr, M, N = cnp.copy2c(arr.T)
    else:
        c_arr, M, N = cnp.copy2c(arr)

    # calculate the lag
    rho = __rho__.lagNRho(ct.byref(c_arr), M, N, ct.c_int(lag))
    cnp.free(c_arr, M, N)

    # return the python-friendly value
    return float(rho)

def RofM(arr, axis:int=1):
    """calculate the Nth lag of 2D input matrix
    
    Parameters:
    ----
    arr: the array over which the Nth lag coherence is being calculated
    lag: the lag index
    axis: the axis corresponding to the elements

    Returns:
    ----
    rho: normalized correlation coefficient at this lag
    """

    # only allow 2D input arrays
    if np.ndim(arr) != 2:
        raise ValueError("Input array must be 2D")
    
    # Make the sample axis the fast axis
    if axis == 1:
        c_arr, M, N = cnp.copy2c(arr.T)
        nele = arr.shape[1]
    else:
        c_arr, M, N = cnp.copy2c(arr)
        nele = arr.shape[0]

    # make a buffer
    rhos = np.zeros(nele-1)

    for lag in range(nele-1):
        # calculate the lag
        rho = __rho__.lagNRho(ct.byref(c_arr), M, N, ct.c_int(lag))
        rhos[lag] = float(rho)
    cnp.free(c_arr, M, N)

    # return the python-friendly value
    return rhos