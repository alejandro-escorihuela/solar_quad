#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 21-11-2021
# alex
# ordre_solar.py

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
    titols = "Sistema Solar"
    # T, N = 1e4, (10**np.linspace(2.0, 6.0, 25)).astype(int)
    # T, N = 1e4, (10**np.linspace(0.1, 4.5, 15)).astype(int)
    T, N = 2e5, (10**np.linspace(0.1, 5.0, 25)).astype(int)
    # quins = [1, 0, 0, 0, 0, 0, 1]
    # quins = [1, 0, 0, 0, 0, 0, 1, 1]
    # quins = [1, 0, 0, 0, 0, 1, 1, 1, 1, 1]
    quins = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    zini, params, planetes = iniSS(quins)
    print("Sistema Solar: " + ", ".join(planetes))

    errkv, navav = [], []
    ri = 0
    n_m = len(cofABA)
    for rk in range(n_m):
        a, b = cofABA[rk]
        x, y = proABA[rk][:len(proABA[rk])//2], proABA[rk][len(proABA[rk])//2:]
        errkv_e, navav_e = [], []
        Nm_ant = -1
        for i in range(len(N)):
            s = len(b)
            n = N[i]
            Nm = int(np.ceil(n/s))
            if Nm > Nm_ant:               
                h = T/Nm
                lerH, lerQ = evolABA_errHQ(zini.copy(), params, Nm, h, a, b, x, y)
                # lerH, lerQ = evolABA_errH(zini.copy(), params, Nm, h, a, b, x, y), np.nan
                errkv_e.append(lerQ)
                # navav_e.append(s*Nm)
                navav_e.append(s/h)
                print("NIA[%02d] = %12.12s -> (ex = %03d/%03d, s = %02d, h = %e, Nm = %07d, Nava = %07d, s/h = %.4e, Tf = %.3e, log10(erH) = %7.3f, log10(erQ) = %7.3f)"
                      % (rk, noms[rk % len(noms)], i + 1, len(N), s, h, Nm, s*Nm, s/h, h*Nm, lerH, lerQ))
                Nm_ant = Nm
        errkv.append(errkv_e)
        navav.append(navav_e)
        ri += 1
        
    ## Gràfiques
    print("Preparant les gràfiques")
    for i in range(len(errkv)):
        plt.plot(np.log10(navav[i]), errkv[i], "C" + str(i%10) + "--", marker = "o", markersize=7.5, linewidth=1.5, label = noms[i % len(noms)])
    plt.title(titols)
    plt.legend()
    plt.savefig("out.pdf", format='pdf')
    plt.show()

