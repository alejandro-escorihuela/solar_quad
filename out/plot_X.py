#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 21-07-2025
# alex
# plot_X.py

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

def linia2list(lin):
    return list(np.array(lin.replace("\n", "").split(" ")).astype(np.float128))

def difnorm2(v1, v2):
    nu, de = 0.0, 0.0
    for i in range(len(v1)):
        nu += (v1[i] - v2[i])*(v1[i] - v2[i])
        de += v2[i]*v2[i];
    return np.sqrt(nu/de)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Hi ha que passar dos fitxers.")
        print("El segon ha de tenir la solució 'exacta'. Les h del primer han d'estar incloses en el segon.")
        exit(-1)
    arxius = sys.argv[1:]
    for i in range(len(arxius)):
        if not os.path.isfile(arxius[i]):
            print("El fitxer", arxius[i], "no existeix.")
            exit(-2)
            
    f1, f2 = open(arxius[0]), open(arxius[1])
    lin1 = list(map(linia2list, f1.readlines()))
    lin2 = list(map(linia2list, f2.readlines()))
    j = 0
    ini = 4 # 10 si no contem el sol
    npl = (len(lin1[0]) - ini)//6
    x, y = [], [[] for _ in range(npl)]
    for i in range(len(lin1)):
        while lin2[j][0] < lin1[i][0]:
            j += 1
        x.append(lin1[i][0])
        for k in range(npl):
            pos = ini + 3*k
            y[k].append(difnorm2(lin1[i][pos:pos + 3], lin2[j][pos:pos + 3]))
    f1.close()
    f2.close()
    
    for i in range(npl):
        plt.plot(np.log10(x), np.log10(y[i]), linewidth = 1.5, label = "Planeta " + str(i))
    plt.legend()
    plt.savefig("erX.pdf", format='pdf')
    plt.show()
        
