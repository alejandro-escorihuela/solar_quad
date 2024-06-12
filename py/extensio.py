#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 05-06-2024
# alex
# extensio.py

def estendreA(cof):
    s = len(cof) + 1
    if s % 2 != 0:
        sep = (s - 1)//2
        a, b = list(cof[:sep]), list(cof[sep:])
        a.append(0.5 - sum(a))
        b.append(1.0 - 2*sum(b))
        a = a + a[::-1]
        b = b + b[-2::-1]
    else:
        sep = s//2
        a, b = list(cof[:sep]), list(cof[sep:])
        a.append(1.0 - 2*sum(a))
        b.append(0.5 - sum(b))
        a = a + a[-2::-1]
        b = b + b[::-1]
    return (a, b)
