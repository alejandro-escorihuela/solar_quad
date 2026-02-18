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
    if len(sys.argv) < 3:
        print("Hi ha que passar més de dos fitxers.")
        print("El darrer ha de tenir la solució 'exacta'. Totes les 'h' han d'estar incloses en l'últim.")
        exit(-1)
    arxius = sys.argv[1:][::-1]
    f = []
    for i in range(len(arxius)):
        if not os.path.isfile(arxius[i]):
            print("El fitxer", arxius[i], "no existeix.")
            exit(-2)
        else:
            f.append(open(arxius[i]))

    lin = []
    for i in range(len(f)):
        lin.append(list(map(linia2list, f[i].readlines())))
        
    x_tot, y_tot = [], []
    for i in range(1, len(lin)):
        x, y = [], []
        j = 0
        for k in range(len(lin[i])):
            while lin[0][j][0] < lin[i][k][0]:
                j += 1
            if abs(lin[0][j - 1][0] - lin[i][k][0]) < 1e-8:
                x.append(lin[0][j - 1][0])
                y.append(difnorm2(lin[0][j - 1][4:6], lin[i][k][4:6]))
                # print(i, lin[0][j - 1][0], lin[i][k][0], lin[0][j - 1][4:6], lin[i][k][4:6])
        x_tot.append(x)
        y_tot.append(y)
    
    for i in range(len(x_tot)):
        plt.plot(np.log10(x_tot[i]), np.log10(y_tot[i]), linewidth = 1.5, label = arxius[i + 1])
    plt.legend()
    plt.savefig("erX.pdf", format='pdf')
    plt.show()
        
