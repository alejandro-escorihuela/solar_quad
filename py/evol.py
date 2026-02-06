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
kpert = ct.CDLL(os.path.dirname(__file__) + "/../c/kpert.so", mode = ct.RTLD_GLOBAL)
evol = ct.CDLL(os.path.dirname(__file__) + "/../c/evol.so", mode = ct.RTLD_GLOBAL)

evol.evolABAsolar.argtypes = (ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_char_p), ct.c_int, ct.c_int, ct.c_longdouble, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.c_char_p, ct.c_int)
evol.evolABAsolar_errH.argtypes = (ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_char_p), ct.c_int, ct.c_int, ct.c_longdouble, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_longdouble))
evol.evolABAsolar_errQ.argtypes = (ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_char_p), ct.c_int, ct.c_int, ct.c_longdouble, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_longdouble))
evol.evolABAsolar_errHQ.argtypes = (ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_char_p), ct.c_int, ct.c_int, ct.c_longdouble, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_longdouble), ct.POINTER(ct.c_longdouble))

evol.evolABAsolar.restype = ct.c_void_p
evol.evolABAsolar_errH.restype = ct.c_void_p
evol.evolABAsolar_errQ.restype = ct.c_void_p
evol.evolABAsolar_errHQ.restype = ct.c_void_p

evol.evolABAkpert.argtypes = (ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_char_p), ct.c_int, ct.c_int, ct.c_longdouble, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.c_char_p, ct.c_int)
evol.evolABAkpert_errH.argtypes = (ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_char_p), ct.c_int, ct.c_int, ct.c_longdouble, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_longdouble))
evol.evolABAkpert_errQ.argtypes = (ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_char_p), ct.c_int, ct.c_int, ct.c_longdouble, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_longdouble))
evol.evolABAkpert_errHQ.argtypes = (ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_char_p), ct.c_int, ct.c_int, ct.c_longdouble, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_char_p), ct.c_int, ct.POINTER(ct.c_longdouble), ct.POINTER(ct.c_longdouble))

evol.evolABAkpert.restype = ct.c_void_p
evol.evolABAkpert_errH.restype = ct.c_void_p
evol.evolABAkpert_errQ.restype = ct.c_void_p
evol.evolABAkpert_errHQ.restype = ct.c_void_p

fun_evo = [evol.evolABAsolar, evol.evolABAkpert]
fun_err = [
    [evol.evolABAsolar_errH, evol.evolABAsolar_errQ, evol.evolABAsolar_errHQ],
    [evol.evolABAkpert_errH, evol.evolABAkpert_errQ, evol.evolABAkpert_errHQ]
]

def evolABA(prob, z, par, Nm, h, a, b, x, y, nom_fit, punts):
    global evol
    np, nz = len(par), len(z)
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
    nom_fit_c = ct.c_char_p(nom_fit.encode())
    fun_evo[prob](z_c, nz, par_c, np, Nm, h, a_c, b_c, len(b), x_c, y_c, len(x), nom_fit_c, punts)
    del par_c, z_c, a_c, b_c, x_c, y_c, nom_fit_c
    
def evolABA_err(prob, err, z, par, Nm, h, a, b, x, y):
    global evol
    np, nz = len(par), len(z)
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
    fun_err[prob][err](z_c, nz, par_c, np, Nm, h, a_c, b_c, len(b), x_c, y_c, len(x), ct.byref(erH), ct.byref(erQ))
    ret_val = [erH.value, erQ.value, (erH.value, erQ.value)]
    del par_c, z_c, a_c, b_c, x_c, y_c
    return ret_val[err]

def evolABAsolar(z, par, Nm, h, a, b, x, y, nom_fit, punts):
    evolABA(0, z, par, Nm, h, a, b, x, y, nom_fit, punts)

def evolABAsolar_errH(z, par, Nm, h, a, b, x, y):
    return evolABA_err(0, 0, z, par, Nm, h, a, b, x, y)

def evolABAsolar_errQ(z, par, Nm, h, a, b, x, y):
    return evolABA_err(0, 1, z, par, Nm, h, a, b, x, y)

def evolABAsolar_errHQ(z, par, Nm, h, a, b, x, y):
    return evolABA_err(0, 2, z, par, Nm, h, a, b, x, y)

def evolABAkpert(z, par, Nm, h, a, b, x, y, nom_fit, punts):
    evolABA(1, z, par, Nm, h, a, b, x, y, nom_fit, punts)

def evolABAkpert_errH(z, par, Nm, h, a, b, x, y):
    return evolABA_err(1, 0, z, par, Nm, h, a, b, x, y)

def evolABAkpert_errQ(z, par, Nm, h, a, b, x, y):
    return evolABA_err(1, 1, z, par, Nm, h, a, b, x, y)

def evolABAkpert_errHQ(z, par, Nm, h, a, b, x, y):
    return evolABA_err(1, 2, z, par, Nm, h, a, b, x, y)
