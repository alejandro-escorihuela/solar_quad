#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 19-07-2025
# alex
# error_solar.py

import matplotlib.pyplot as plt
from matplotlib import rc
# from decimal import *
import numpy as np
import sys
sys.path.insert(0, './py')
from metsNIA import *
from kpert import *
from evol import *

# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size' : 18})
# rc('text', usetex=True)
# plt.rcParams['figure.dpi'] = 288

if __name__ == "__main__":
    T, h, punts = 500*2*np.pi, 2.5e-2, int(1e4)
    Nfixes = 3*T/h
    exc, eps, alp = 0.8, 1e-4, 1.0
    zini, params = inikpert(exc, eps, alp)
    zini = ["0.20000262922934886", "-2.1607753661284597e-06", "8.048436279808638e-06", "2.999880962197589"]
    titols = "Kepler pertorbat ϵ = %s, α = %s, e = %f" % (params[0], params[1], exc)
    print(titols)
    print("z₀ =", zini)    
    for rk in range(len(cofABA)):
        a, b = cofABA[rk]
        x, y = proABA[rk][:len(proABA[rk])//2], proABA[rk][len(proABA[rk])//2:]
        s = len(b)
        h = s*T/Nfixes
        Nm = int(T//h)
        punts = Nm
        nom_fit = "./out/" + noms[rk % len(noms)] + ".txt"
        evolABAkpert(zini.copy(), params, Nm, float(h), a, b, x, y, nom_fit, punts)
        print("NIA[%02d] -> (arxiu = %s, punts = %d, h = %e, Nm = %d, Nava = %d, s/h = %.4e, Tf = %.3e)"
                      % (rk, nom_fit, punts, h, Nm, s*Nm, s/h, T))

