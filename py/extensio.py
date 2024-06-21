#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 05-06-2024
# alex
# extensio.py

import numpy as np

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
