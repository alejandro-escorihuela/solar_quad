#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 21-11-2021
# alex
# ordre_solar.py

import matplotlib.pyplot as plt
from matplotlib import rc
#from decimal import *
import multiprocessing as mupr
import numpy as np
import sys
sys.path.insert(0, './py')
from metsNIA import *
from solar import *
from evol import *

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size' : 18})
#rc('text', usetex=True)
#plt.rcParams['figure.dpi'] = 288

def prc_evol(val, pid, valors):
    zini, rk, Nm, h, nex, hoq, par = val
    a, b = cofABA[rk]
    x, y = proABA[rk][:len(proABA[rk])//2], proABA[rk][len(proABA[rk])//2:]
    if hoq == 0:
        lerH, lerQ = evolABAsolar_errH(zini, par, Nm, h, a, b, x, y), np.nan
    elif hoq == 1:
        lerH, lerQ = evolABAsolar_errHQ(zini, par, Nm, h, a, b, x, y)
    valors[pid] = [rk, nex, len(b), Nm, h, pid, lerH, lerQ]

def ordenar_par(e):
    return e[1]*e[2] # Nava = s*Nm

def ordenar_res(e):
    return e[0]

def res1a2(lis):
    nav, err = [[] for _ in range(len(lis))], [[] for _ in range(len(lis))]
    for i in range(len(lis)):
        for j in range(len(lis[i])):
            nav[i].append(lis[i][j][0])
            err[i].append(lis[i][j][1])
    return nav, err

def print_info(pr, rk, nex, n, Nm, s, h, err1, err2):
    print("pr = %03d | NIA[%02d] = %12.12s -> (ex = %03d/%03d, s = %02d, h = %e, Nm = %07d, Nava = %07d, s/h = %.4e, Tf = %.3e, log10(erH) = %7.3f, log10(erQ) = %7.3f)"
          % (pr, rk, noms[rk % len(noms)], nex, n, s, h, Nm, s*Nm, s/h, h*Nm, err1, err2)) 

if __name__ == "__main__":
    max_prcs = 25
    manager = mupr.Manager()
    valpr = manager.dict()      
    T, N = 1e4, (10**np.linspace(2.0, 6.0, 25)).astype(int)
    # T, N = 1e4, (10**np.linspace(0.1, 4.5, 15)).astype(int)
    # T, N = 2e5, (10**np.linspace(0.1, 5.0, 25)).astype(int)
    # quins = [1, 0, 0, 0, 0, 0, 1]
    # quins = [1, 0, 0, 0, 0, 0, 1, 1]
    # quins = [1, 0, 0, 0, 0, 1, 1, 1, 1, 1]
    quins = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    zini, params, planetes = iniSS(quins)
    hoq = 1
    titols = "Sistema Solar"
    print(titols + ": " + ", ".join(planetes))

    # Preparem un vector de paràmetres per a executar en paraŀlel
    v_par = []
    ri = 0
    n_m = len(cofABA)
    for rk in range(n_m):
        Nm_ant = -1
        for i in range(len(N)):
            s = len(cofABA[rk][1])
            n = N[i]
            Nm = int(np.ceil(n/s))
            if Nm > Nm_ant:               
                h = T/Nm
                v_par.append([rk, s, Nm, h, i + 1])
                Nm_ant = Nm
    v_par.sort(key = ordenar_par)

    # Executem en paraŀlel
    lis_res = [[] for _ in range(n_m)]
    print("Nombre d'execucions total:      ", len(v_par))
    print("Nombre d'execucions en paraŀlel:", max_prcs)
    for i in range(len(v_par)//max_prcs + 1):
        prcs = []
        valpr.clear() # Necessari si pr no té sempre els mateixos valors
        for pr in range(i*max_prcs, min((i + 1)*max_prcs, len(v_par))):
            x = mupr.Process(target = prc_evol, args=([zini.copy(), v_par[pr][0], v_par[pr][2], v_par[pr][3], v_par[pr][4], hoq, params], pr, valpr))
            prcs.append(x)
            x.daemon = True
            x.start()
        for prc in prcs:
            prc.join()
        for clau in valpr:
            rk, nex, s, Nm, h, pr, lerH, lerQ = valpr[clau]
            if hoq == 0:
                lis_res[rk].append([s/h, lerH])
            elif hoq == 1:
                lis_res[rk].append([s/h, lerQ])
            print_info(pr, rk, nex, len(N), Nm, s, h, lerH, lerQ)    

    # Organitzem les dades per a les gràfiques
    for i in range(len(lis_res)):
       lis_res[i].sort(key = ordenar_res)
    navav, errkv = res1a2(lis_res)


    # Dibuixem les gràfiques
    print("Preparant les gràfiques")
    if (len(N) <= 25):
        for i in range(len(errkv)):
            plt.plot(np.log10(navav[i]), errkv[i], "C" + str(i%10) + "--", marker = "o", markersize=7.5, linewidth=1.5, label = noms[i % len(noms)])
    else:
        for i in range(len(errkv)):
            plt.plot(np.log10(navav[i]), errkv[i], "C" + str(i%10), label = noms[i % len(noms)])             
    plt.title(titols)
    plt.xlabel(r"$\log_{10}(s/h)$")
    if hoq == 0:
        plt.ylabel(r"$\log_{10}$" + "(error energies)")
    elif hoq == 1:
        plt.ylabel(r"$\log_{10}$" + "(error posicions)")
    plt.legend()
    plt.savefig("out_par.pdf", format='pdf')
    plt.show()

