#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 05-06-2024
# alex
# evol.py

import numpy as np
import ctypes as ct
import os
from solar import *

vector = ct.CDLL(os.path.dirname(__file__) + "/../c/vector.so", mode = ct.RTLD_GLOBAL)
solar = ct.CDLL(os.path.dirname(__file__) + "/../c/solar.so", mode = ct.RTLD_GLOBAL)
evol = ct.CDLL(os.path.dirname(__file__) + "/../c/evol.so", mode = ct.RTLD_GLOBAL)
evol.evolABAsolar_errH.argtypes = (ct.POINTER(ct.c_longdouble), ct.POINTER(ct.c_longdouble), ct.c_int, ct.c_int, ct.c_longdouble, ct.POINTER(ct.c_longdouble), ct.POINTER(ct.c_longdouble), ct.c_int, ct.POINTER(ct.c_longdouble), ct.POINTER(ct.c_longdouble), ct.c_int, ct.POINTER(ct.c_longdouble))    
evol.evolABAsolar_errQ.argtypes = (ct.POINTER(ct.c_longdouble), ct.POINTER(ct.c_longdouble), ct.c_int, ct.c_int, ct.c_longdouble, ct.POINTER(ct.c_longdouble), ct.POINTER(ct.c_longdouble), ct.c_int, ct.POINTER(ct.c_longdouble), ct.POINTER(ct.c_longdouble), ct.c_int, ct.POINTER(ct.c_longdouble))    
evol.evolABAsolar_errHQ.argtypes = (ct.POINTER(ct.c_longdouble), ct.POINTER(ct.c_longdouble), ct.c_int, ct.c_int, ct.c_longdouble, ct.POINTER(ct.c_longdouble), ct.POINTER(ct.c_longdouble), ct.c_int, ct.POINTER(ct.c_longdouble), ct.POINTER(ct.c_longdouble), ct.c_int, ct.POINTER(ct.c_longdouble), ct.POINTER(ct.c_longdouble))    
evol.evolABAsolar_errH.restype = ct.c_void_p
evol.evolABAsolar_errQ.restype = ct.c_void_p
evol.evolABAsolar_errHQ.restype = ct.c_void_p

def evolABA_errH(z, Nm, h, a, b, x, y):
    masses_c = ss_ma.ctypes.data_as(ct.POINTER(ct.c_longdouble))
    z_c = z.ctypes.data_as(ct.POINTER(ct.c_longdouble))
    a_c = np.array(a).astype("float128").ctypes.data_as(ct.POINTER(ct.c_longdouble))
    b_c = np.array(b).astype("float128").ctypes.data_as(ct.POINTER(ct.c_longdouble))
    x_c = np.array(x).astype("float128").ctypes.data_as(ct.POINTER(ct.c_longdouble))
    y_c = np.array(y).astype("float128").ctypes.data_as(ct.POINTER(ct.c_longdouble))
    erH = ct.c_longdouble(0.0)
    evol.evolABAsolar_errH(z_c, masses_c, ss_np, Nm, h, a_c, b_c, len(b), x_c, y_c, len(x), ct.byref(erH))
    return erH.value

def evolABA_errQ(z, Nm, h, a, b, x, y):
    masses_c = ss_ma.ctypes.data_as(ct.POINTER(ct.c_longdouble))
    z_c = z.ctypes.data_as(ct.POINTER(ct.c_longdouble))
    a_c = np.array(a).astype("float128").ctypes.data_as(ct.POINTER(ct.c_longdouble))
    b_c = np.array(b).astype("float128").ctypes.data_as(ct.POINTER(ct.c_longdouble))
    x_c = np.array(x).astype("float128").ctypes.data_as(ct.POINTER(ct.c_longdouble))
    y_c = np.array(y).astype("float128").ctypes.data_as(ct.POINTER(ct.c_longdouble))
    erQ = ct.c_longdouble(0.0)
    evol.evolABAsolar_errQ(z_c, masses_c, ss_np, Nm, h, a_c, b_c, len(b), x_c, y_c, len(x), ct.byref(erQ))
    return erQ.value

def evolABA_errHQ(z, Nm, h, a, b, x, y):
    masses_c = ss_ma.ctypes.data_as(ct.POINTER(ct.c_longdouble))
    z_c = z.ctypes.data_as(ct.POINTER(ct.c_longdouble))
    a_c = np.array(a).astype("float128").ctypes.data_as(ct.POINTER(ct.c_longdouble))
    b_c = np.array(b).astype("float128").ctypes.data_as(ct.POINTER(ct.c_longdouble))
    x_c = np.array(x).astype("float128").ctypes.data_as(ct.POINTER(ct.c_longdouble))
    y_c = np.array(y).astype("float128").ctypes.data_as(ct.POINTER(ct.c_longdouble))
    erH = ct.c_longdouble(0.0)
    erQ = ct.c_longdouble(0.0)
    evol.evolABAsolar_errHQ(z_c, masses_c, ss_np, Nm, h, a_c, b_c, len(b), x_c, y_c, len(x), ct.byref(erH), ct.byref(erQ))
    return (erH.value, erQ.value)
