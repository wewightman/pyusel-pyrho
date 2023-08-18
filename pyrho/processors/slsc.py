"""SLSC processing wrapper for a given acq"""
import ctypes as ct
from cinpy import copy2c, free
import numpy as np
from pyrho.trig import c_norm_engine, c_pw_engine

from abc import ABC, abstractmethod

class TXType(ABC):
    @classmethod
    @abstractmethod
    def gentabs(self, points):
        raise NotImplementedError
    
    @classmethod
    @abstractmethod
    def c_gentabs(self, points, Np:int):
        raise NotImplementedError
    
    @classmethod
    @abstractmethod
    def cleartabs(self):
        raise NotImplementedError
    
class RXType(ABC):
    @classmethod
    @abstractmethod
    def gentabs(self, points):
        raise NotImplementedError
    
    @classmethod
    @abstractmethod
    def c_gentabs(self, points, Np:int):
        raise NotImplementedError
    
    @classmethod
    @abstractmethod
    def cleartabs(self):
        raise NotImplementedError
    
class FullRX(RXType):
    def __init__(self, xrefs, c=1540, dtype=ct.c_float):
        # vaidate the inputs are valid shapes
        if np.ndim(xrefs) != 2: raise ValueError("xrefs must be 2D")
        if xrefs.shape[1] != 3: raise ValueError("xrefs must be Nacq by 3")
        self.Nrx = xrefs.shape[0]

        # copy attributes and convert data to lists of ctypes
        self.xrefs = [ct.cast((dtype * 3)(*[dtype(x) for x in xrefs[irx][:]]), ct.POINTER(dtype))
                      for irx in range(self.Nrx)]
        self.c = c
        self.dtype = dtype

        # initialize and fill RXtabs
        self.RXtabs = None

    def gentabs(self, points):
        if self.RXtabs is not None:
            self.cleartabs()
        
        if np.ndim(points) != 2: raise ValueError("points matrix must be 2D")
        if points.shape[1] != 3: raise ValueError("points matrix must be P by 3")

        self.Np = points.shape[0]

        c_points, Mp, Np = copy2c(points)

        self.RXtabs = []
        for itx in range(self.Nrx):
            tabs = c_norm_engine(
                self.xrefs[itx],
                c_points,
                Mp.value,
                self.c,
                self.dtype
            )
            self.RXtabs.append(tabs)

        free(c_points, Mp, Np)

    def c_gentabs(self, c_points, Np:int):
        if self.RXtabs is not None:
            self.cleartabs()

        self.RXtabs = []
        for itx in range(self.Nrx):
            tabs = c_norm_engine(
                self.xrefs[itx],
                c_points,
                Np,
                self.c,
                self.dtype
            )
            self.RXtabs.append(tabs)
    
    def cleartabs(self):
        if self.RXtabs is not None:
            for tabs in self.RXtabs:
                free(tabs, ct.c_int(self.Np))
            self.TXtabs = None
            self.Np = None

    def __del__(self):
        # clear all tabs if any exists
        self.cleartabs()

class PWTX(TXType):
    def __init__(self, alphas, xrefs, trefs, betas=None, c=1540, dtype=ct.c_float):
        # vaidate the inputs are valid shapes
        if np.ndim(xrefs) != 2: raise ValueError("xrefs must be 2D")
        if xrefs.shape[1] != 3: raise ValueError("xrefs must be Nacq by 3")
        self.Ntx = xrefs.shape[0]
        if betas is None: betas = np.zeros(self.Ntx)
        if np.ndim(alphas) != 1: raise ValueError("alphas must be a 1D vector")
        if np.ndim(betas) != 1: raise ValueError("betas must be a 1D vector")
        if np.ndim(trefs) != 1: raise ValueError("trefs must be a 1D vector")
        if len(alphas) != self.Ntx: raise ValueError("alphas must have the same length as xrefs")
        if len(betas) != self.Ntx: raise ValueError("betas must have the same length as xrefs")
        if len(trefs) != self.Ntx: raise ValueError("trefs must have the same length as xrefs")

        # initialize and fill TXtabs
        self.TXtabs = None
        self.Np = None

        # copy attributes and convert data to lists of ctypes
        self.norms = [ct.cast((dtype * 3)(
                        np.sin(alpha), 
                        np.cos(alpha)*np.cos(beta), 
                        np.sin(beta)), ct.POINTER(dtype))
                      for alpha, beta in zip(alphas, betas)]
        self.xrefs = [ct.cast((dtype * 3)(*[dtype(x) for x in xrefs[itx][:]]), ct.POINTER(dtype))
                      for itx in range(self.Ntx)]
        self.trefs = [dtype(tref) for tref in trefs]
        self.c = c
        self.dtype = dtype

    def gentabs(self, points):
        if self.TXtabs is not None:
            self.cleartabs()
        
        if np.ndim(points) != 2: raise ValueError("points matrix must be 2D")
        if points.shape[1] != 3: raise ValueError("points matrix must be P by 3")
        self.Np = points.shape[0]
        c_points, Mp, Np = copy2c(points)

        self.TXtabs = []
        for itx in range(self.Ntx):
            tabs = c_pw_engine(
                self.xrefs[itx],
                self.trefs[itx],
                self.norms[itx],
                c_points,
                Mp.value,
                self.c,
                self.dtype
            )
            self.TXtabs.append(tabs)

        free(c_points, Mp, Np)

    def c_gentabs(self, c_points, Np:int):
        if self.TXtabs is not None:
            self.cleartabs()

        self.TXtabs = []
        for itx in range(self.Ntx):
            tabs = c_pw_engine(
                self.xrefs[itx],
                self.trefs[itx],
                self.norms[itx],
                c_points,
                Np,
                self.c,
                self.dtype
            )
            self.TXtabs.append(tabs)
    
    def cleartabs(self):
        if self.TXtabs is not None:
            for tabs in self.TXtabs:
                free(tabs, ct.c_int(self.Np))
            self.TXtabs = None
            self.Np = None

    def __del__(self):
        # clear all tabs if any exists
        self.cleartabs()

class SLSCProc():
    def __init__(self, c, points, 
                 tx:TXType,
                 rx:RXType,
                 lags:list|int=1, 
                 fnum:list|float|None=None,
                 fnorm:list|np.ndarray|None=None,
                 dtype=ct.c_float, **kwargs):
        """Initialize a SLSC processor. 
        Define the speed of sound, the 3D points to be reconstructed, the effective fnumber(s) and relative axis(axes)

        Parameters:
        ----
        c: the assumed speed of sound in the medium in m/s
        tx: the transmit object
        points: a P by 3 matrix where each row represents the spatial coordinates of a reconstruction point in m
        lags: list of lag indices to integrate over
        fnum: fnumber(s) to use when reconstructucting effective apertures. If None (default), will use all channels for all points.
        fnorm: normal vectors
        """

        # check that the correct transmission parameters are included
        self.tx = tx
        self.rx = rx
        self.c = dtype(c)
        self.points, self.c_Npoints, _ = copy2c(points)
        self.lags = lags
        self.fnum = fnum
        self.fnorm = fnorm
        self.dtype = dtype

        self.tx.c_gentabs(self.points, self.c_Npoints.value)
        self.rx.c_gentabs(self.points, self.c_Npoints.value)

    def __call__(self, data):
        print("oopa")
