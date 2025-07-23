#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 21-07-2025
# alex
# plot_H.py

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

def linia2list(lin):
    return list(np.array(lin.replace("\n", "").split(" ")).astype(np.float128))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Hi ha que passar com a mínim un fitxer.")
        exit(-1)
    arxius = sys.argv[1:]
    for i in range(len(arxius)):
        if not os.path.isfile(arxius[i]):
            print("El fitxer", arxius[i], "no existeix.")
            exit(-2)
            
    x, y = [], []
    for i in range(len(arxius)):
        x_i, y_i = [], []
        with open(arxius[i]) as f:
            linies = f.readlines()    
            for j in range(len(linies)):
                nums = linia2list(linies[j])
                x_i.append(nums[0])
                y_i.append(nums[3])
        x.append(x_i)
        y.append(y_i)

    for i in range(len(arxius)):
        plt.plot(np.log10(x[i]), np.log10(y[i]), linewidth = 1.5, label = arxius[i])
    plt.legend()
    plt.savefig("erH.pdf", format='pdf')
    plt.show()
