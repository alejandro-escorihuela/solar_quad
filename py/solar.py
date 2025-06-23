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
    mas = list(map(str, m_sel))
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

