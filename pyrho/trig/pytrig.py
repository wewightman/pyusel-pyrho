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

