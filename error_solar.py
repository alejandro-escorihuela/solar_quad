#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 19-07-2025
# alex
# error_solar.py

import matplotlib.pyplot as plt
from matplotlib import rc
#from decimal import *
import numpy as np
import sys
sys.path.insert(0, './py')
from metsNIA import *
from solar import *
from evol import *

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size' : 18})
#rc('text', usetex=True)
#plt.rcParams['figure.dpi'] = 288

if __name__ == "__main__":
    T, h, punts = 2e5, 87.95/2, 2000
    T, h, punts = 2e6, 20, 10000
    # T, h, punts = 2e6, 6, int(2e6//6)
    Nfixes = 3*T/h
    # quins = [1, 0, 0, 0, 0, 0, 1]
    # quins = [1, 0, 0, 0, 0, 0, 1, 1]
    # quins = [1, 0, 0, 0, 0, 1, 1, 1, 1, 1]
    quins = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    zini, params, planetes = iniSS(quins)
    print("Sistema Solar: " + ", ".join(planetes))
    for rk in range(len(cofABA)):
        a, b = cofABA[rk]
        x, y = proABA[rk][:len(proABA[rk])//2], proABA[rk][len(proABA[rk])//2:]
        s = len(b)
        h = s*T/Nfixes
        Nm = int(T//h)
        nom_fit = "./out/" + noms[rk % len(noms)] + ".txt"
        evolABAsolar(zini.copy(), params, Nm, float(h), a, b, x, y, nom_fit, punts)
        print("NIA[%02d] -> (arxiu = %s, punts = %d, h = %e, Nm = %d, Nava = %d, s/h = %.4e, Tf = %.3e)"
                      % (rk, nom_fit, punts, h, Nm, s*Nm, s/h, T))

