#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 06-02-2026
# alex
# kpert.py

import numpy as np
import mpmath as mp

def inikpert(e, eps, alp):
    mp.mp.dps = 50
    mp_e = mp.mpf(str(e))
    mp_z = [mp.mpf("1") - mp_e, mp.mpf("0"), mp.mpf("0"), mp.sqrt((mp.mpf("1") + mp_e)/(mp.mpf("1") - mp_e))]
    str_z = [mp.nstr(zi, 50) for zi in mp_z]
    return str_z, [str(eps), str(alp)]
