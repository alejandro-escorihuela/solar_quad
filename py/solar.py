#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 05-06-2024
# alex
# solar.py

import numpy as np

# Masses sistema solar
ss_ma = np.array([1.0,
                  1.660130719E-07, 2.447712418E-06, 3.002609351E-06, 3.226998492E-07,
                  9.543137255E-04, 2.857310206E-04,
                  4.364519859E-05, 5.148818502E-05, 6.571141277E-09]).astype("float128")
#ss_ma = list(map(str, ss_ma))

# Nombre de planetes sistema solar
ss_np = len(ss_ma)

def iniSS():
    q, p = np.zeros((ss_np, 3)), np.zeros((ss_np, 3))
    ic = 0
    # Sol
    q[ic][0] = 0.0
    q[ic][1] = 0.0
    q[ic][2] = 0.0
    p[ic][0] = 0.0
    p[ic][1] = 0.0
    p[ic][2] = 0.0
    # Mercuri
    ic += 1
    q[ic][0] = -3.859054833019958E-01
    q[ic][1] = -1.654307689835696E-03
    q[ic][2] = 3.481192038201179E-02
    p[ic][0] = np.float128(ss_ma[ic])*-5.294168138166344E-03
    p[ic][1] = np.float128(ss_ma[ic])*-2.691453995673129E-02
    p[ic][2] = np.float128(ss_ma[ic])*-1.714386730721826E-03
    # Venus 
    ic += 1
    q[ic][0] = 7.293167005709342E-02
    q[ic][1] = -7.175590811561052E-01
    q[ic][2] = -1.415191215112003E-02
    p[ic][0] = np.float128(ss_ma[ic])*1.998700269409145E-02
    p[ic][1] = np.float128(ss_ma[ic])*1.911135291296663E-03
    p[ic][2] = np.float128(ss_ma[ic])*-1.127428703200971E-03
    # La Terra
    ic += 1
    q[ic][0] = -1.734173457390217E-01
    q[ic][1] = 9.736937590796603E-01
    q[ic][2] = -1.582107821713564E-04
    p[ic][0] = np.float128(ss_ma[ic])*-1.720709737865684E-02
    p[ic][1] = np.float128(ss_ma[ic])*-3.125258586509626E-03
    p[ic][2] = np.float128(ss_ma[ic])*-1.120932427483096E-07
    # Mart
    ic += 1
    q[ic][0] = -1.581869997754351E+00
    q[ic][1] = -3.828978953153063E-01
    q[ic][2] = 3.059389984291728E-02
    p[ic][0] = np.float128(ss_ma[ic])*3.855008143708335E-03
    p[ic][1] = np.float128(ss_ma[ic])*-1.238923952632399E-02
    p[ic][2] = np.float128(ss_ma[ic])*-3.543305537962489E-04
    # Jupiter
    ic += 1
    q[ic][0] = -4.259467773894452E+00
    q[ic][1] = -3.361194945480983E+00
    q[ic][2] = 1.092145047021746E-01
    p[ic][0] = np.float128(ss_ma[ic])*4.586299412789570E-03
    p[ic][1] = np.float128(ss_ma[ic])*-5.564874896880609E-03
    p[ic][2] = np.float128(ss_ma[ic])*-7.945799167252124E-05
    # Saturn
    ic += 1
    q[ic][0] = 4.789734270644876E-02
    q[ic][1] = -1.005701578869786E+01
    q[ic][2] = 1.729539827294794E-01
    p[ic][0] = np.float128(ss_ma[ic])*5.271615539820981E-03
    p[ic][1] = np.float128(ss_ma[ic])*8.862372960977510E-06
    p[ic][2] = np.float128(ss_ma[ic])*-2.100595394827879E-04
    # Urà
    ic += 1
    q[ic][0] = 1.772328745814774E+01
    q[ic][1] = 9.063002917185520E+00
    q[ic][2] = -1.959478058581542E-01
    p[ic][0] = np.float128(ss_ma[ic])*-1.819603624325987E-03
    p[ic][1] = np.float128(ss_ma[ic])*3.318475309448707E-03
    p[ic][2] = np.float128(ss_ma[ic])*3.577108114482214E-05
    # Neptú
    ic += 1
    q[ic][0] = 2.868162693362844E+01
    q[ic][1] = -8.591658348777845E+00
    q[ic][2] = -4.840680053568654E-01
    p[ic][0] = np.float128(ss_ma[ic])*8.802822510921428E-04
    p[ic][1] = np.float128(ss_ma[ic])*3.025692572392946E-03
    p[ic][2] = np.float128(ss_ma[ic])*-8.295671458148408E-05
    # Plutó
    ic += 1
    q[ic][0] = 1.077826511187572E+01
    q[ic][1] = -3.168642408143715E+01
    q[ic][2] = 2.729178542838963E-01
    p[ic][0] = np.float128(ss_ma[ic])*3.030812460422457E-03
    p[ic][1] = np.float128(ss_ma[ic])*3.426619083057393E-04
    p[ic][2] = np.float128(ss_ma[ic])*-9.199095031922107E-04
    z = np.concatenate((q.reshape(ss_np*3), p.reshape(ss_np*3))).astype("float128")
    z = list(map(str, z))
    return z

def expand_masses(vec):
    np = len(vec)
    m, M, mb, nu, muG = list(vec), [0.0]*np, [0.0]*np, [0.0]*np, [0.0]*np
    for i in range(np):
        Mi = 0.0
        for j in range(i + 1):
            Mi += m[j]
        M[i] = Mi
    mb[0] = M[np - 1]
    nu[0] = m[0]/M[0]
    muG[0] = 0.0 # no té validesa física
    for i in range(1, np):
        mb[i] = M[i - 1]*m[i]/M[i]
        nu[i] = m[i]/M[i]
        muG[i] = m[0]*m[i]/mb[i]
    return list(map(str, m + M + mb + nu + muG))
