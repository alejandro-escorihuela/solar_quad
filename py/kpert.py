#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 06-02-2026
# alex
# kpert.py

import numpy as np

def inikpert(e):
    return [1 - e, 0.0, 0.0, np.sqrt((1 + e)/(1 - e))]

