#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 05-06-2024
# alex
# solar.py

import numpy as np

def iniSS(quins):
    ## Exemples
    ## Júpiter - Plutó
    ##   quins = [1, 0, 0, 0, 0, 1, 1, 1, 1, 1]
    ## Mercuri - Mart (tam = 5)
    ##   quins = [1, 1, 1, 1, 1]
    ## Mercuri - Plutó (tam = 10)
    ##   quins = [1]*10
    ## Mercuri - Plutó + satèl·lits + planetes nans (tam = 18)
    ##   quins = [1]*18
    
    ## dades a 28/06/1969 00:00
    m, q, v, noms = llegirfit("./ini/horizons1969.txt")
    
    tam_sel = quins.count(1)
    primerp = quins[1:].index(not 0)
    if quins[0] != 0 and primerp > 0:
        m[0] = sum(m[:primerp + 1])
    q_sel, v_sel = np.zeros((tam_sel, 3)), np.zeros((tam_sel, 3))
    m_sel, noms_sel = [0.0]*tam_sel, [""]*tam_sel
    j = 0
    for i in range(len(quins)):
        if quins[i] != 0:
            m_sel[j], noms_sel[j] = m[i], noms[i]
            q_sel[j][0], q_sel[j][1], q_sel[j][2] = q[i][0], q[i][1], q[i][2]
            v_sel[j][0], v_sel[j][1], v_sel[j][2] = v[i][0], v[i][1], v[i][2]
            j += 1
    z = np.concatenate((q_sel.reshape(tam_sel*3), v_sel.reshape(tam_sel*3))).astype("float128")
    z = list(map(str, z))
    mas = expand_masses(m_sel)
    return z, mas, noms_sel

def llegirfit(fit):
    with open(fit) as f:
        linies = f.readlines()
        tam = len(linies)
        npl = tam - 2
        q, v = np.zeros((npl, 3)), np.zeros((npl, 3))
        m, noms  = [0.0]*npl, [""]*npl
        for i in range(npl):
            lin = np.array(linies[i + 2].replace("\n", "").split(" "))[-7:].astype(np.float128)
            noms[i] = " ".join(linies[i + 2].replace("\n", "").split(" ")[:-7])
            m[i] = lin[0]
            q[i][0], q[i][1], q[i][2] = lin[1], lin[2], lin[3]
            v[i][0], v[i][1], v[i][2] = lin[4], lin[5], lin[6]
    return m, q, v, noms
    
def expand_masses(vec):
    npl = len(vec)
    Gm, GM, Gmb, nu, mu = list(vec).copy(), [0.0]*npl, [0.0]*npl, [0.0]*npl, [0.0]*npl
    GM[0] = Gm[0]
    for i in range(npl):
        GM[i] = GM[i - 1] + Gm[i]
    Gmb[0] = GM[npl - 1]
    nu[0] = Gm[0]/GM[0]
    mu[0] = 0.0
    for i in range(1, npl):
        Gmb[i] = GM[i - 1]*Gm[i]/GM[i]
        nu[i] = Gm[i]/GM[i]
        mu[i] = Gm[0]*Gm[i]/Gmb[i]
    return list(map(str, Gm + GM + Gmb + nu + mu))
