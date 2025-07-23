#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 05-06-2024
# alex
# extensio.py

import numpy as np
import mpmath as mp

def estendreA(cof):
    s = len(cof) - 1
    if s % 2 != 0:
        sep = (s - 1)//2 + 1
        a, b = list(cof[:sep]), list(cof[sep:])
        a = a + a[::-1]
        b = b + b[-2::-1]
    else:
        sep = s//2 + 1
        a, b = list(cof[:sep]), list(cof[sep:])
        a = a + a[-2::-1]
        b = b + b[::-1]
    # a = list(np.array(a, dtype=np.float128))
    # b = list(np.array(b, dtype=np.float128))
    return (a, b)

def str2split(cof):
    mp.mp.dps = 50
    mcof = list(map(mp.mpf, cof))
    alp = []
    for c in mcof:
        alp.append(c/mp.mpf("2"))
        alp.append(c/mp.mpf("2"))
    a, b = [alp[0]], [alp[0] + alp[1]]
    for i in range(1, len(alp)//2):
        a.append(alp[2*i - 1] + alp[2*i])
        b.append(alp[2*i] + alp[2*i + 1])
    a.append(alp[-1])
    astr = [mp.nstr(ai, 50) for ai in a]
    bstr = [mp.nstr(bi, 50) for bi in b]
    return (astr, bstr)    

