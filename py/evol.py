#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 05-06-2024
# alex
# evol.py

import numpy as np
import ctypes as ct
import os

vector = ct.CDLL(os.path.dirname(__file__) + "/../c/vector.so", mode = ct.RTLD_GLOBAL)
solar = ct.CDLL(os.path.dirname(__file__) + "/../c/solar.so", mode = ct.RTLD_GLOBAL)
evol = ct.CDLL(os.path.dirname(__file__) + "/../c/evol.so", mode = ct.RTLD_GLOBAL)

evol.evolABAsolar_errH.argtypes = (ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.c_int, ct.c_longdouble, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_longdouble))
evol.evolABAsolar_errQ.argtypes = (ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.c_int, ct.c_longdouble, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_longdouble))
evol.evolABAsolar_errHQ.argtypes = (ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.c_int, ct.c_longdouble, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_longdouble), ct.POINTER(ct.c_longdouble))

evol.evolABAsolar_errH.restype = ct.c_void_p
evol.evolABAsolar_errQ.restype = ct.c_void_p
evol.evolABAsolar_errHQ.restype = ct.c_void_p

def evolABA_errH(z, par, Nm, h, a, b, x, y):
    global evol
    np = len(par)
    par_c = (ct.c_char_p*np)()
    z_c = (ct.c_char_p*len(z))()
    a_c = (ct.c_char_p*len(a))()
    b_c = (ct.c_char_p*len(b))()
    x_c = (ct.c_char_p*len(x))()
    y_c = (ct.c_char_p*len(y))()
    par_c[:] = [item.encode() for item in par]
    z_c[:] = [item.encode() for item in z]
    a_c[:] = [item.encode() for item in a]
    b_c[:] = [item.encode() for item in b]
    x_c[:] = [item.encode() for item in x]
    y_c[:] = [item.encode() for item in y]
    erH = ct.c_longdouble(0.0)
    evol.evolABAsolar_errH(z_c, par_c, np, Nm, h, a_c, b_c, len(b), x_c, y_c, len(x), ct.byref(erH))
    del par_c, z_c, a_c, b_c, x_c, y_c
    return erH.value

def evolABA_errQ(z, par, Nm, h, a, b, x, y):
    global evol
    np = len(par)
    par_c = (ct.c_char_p*np)()
    z_c = (ct.c_char_p*len(z))()
    a_c = (ct.c_char_p*len(a))()
    b_c = (ct.c_char_p*len(b))()
    x_c = (ct.c_char_p*len(x))()
    y_c = (ct.c_char_p*len(y))()
    par_c[:] = [item.encode() for item in par]
    z_c[:] = [item.encode() for item in z]
    a_c[:] = [item.encode() for item in a]
    b_c[:] = [item.encode() for item in b]
    x_c[:] = [item.encode() for item in x]
    y_c[:] = [item.encode() for item in y]    
    erQ = ct.c_longdouble(0.0)
    evol.evolABAsolar_errQ(z_c, par_c, np, Nm, h, a_c, b_c, len(b), x_c, y_c, len(x), ct.byref(erQ))
    return erQ.value

def evolABA_errHQ(z, par, Nm, h, a, b, x, y):
    global evol
    np = len(par)
    par_c = (ct.c_char_p*np)()
    z_c = (ct.c_char_p*len(z))()
    a_c = (ct.c_char_p*len(a))()
    b_c = (ct.c_char_p*len(b))()
    x_c = (ct.c_char_p*len(x))()
    y_c = (ct.c_char_p*len(y))()
    par_c[:] = [item.encode() for item in par]
    z_c[:] = [item.encode() for item in z]
    a_c[:] = [item.encode() for item in a]
    b_c[:] = [item.encode() for item in b]
    x_c[:] = [item.encode() for item in x]
    y_c[:] = [item.encode() for item in y]
    erH = ct.c_longdouble(0.0)
    erQ = ct.c_longdouble(0.0)
    evol.evolABAsolar_errHQ(z_c, par_c, np, Nm, h, a_c, b_c, len(b), x_c, y_c, len(x), ct.byref(erH), ct.byref(erQ))
    return (erH.value, erQ.value)
