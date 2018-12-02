import numpy as np, scipy as sc, matplotlib as mpl, sys
import matplotlib.pyplot as plt

N = 4
rho = 1000.0
gamma = 0.143
u = 1.0
Length = 1.0
Breadth = 0.5

del_x = Length/N
del_y = Breadth/N

T_e0 = 333.0
T_w0 = 293.0
T_n0 = 303.0
T_s0 = 353.0

A = np.zeros(N**4).reshape(N**2,N**2)
T = np.zeros(N**2)
b = np.zeros(N**2)

n = 1
print l
     
