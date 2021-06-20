# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 09:48:51 2021

@author: Olaf

Chapter 2
"""

import numpy as np
from random import random
import matplotlib.pyplot as plt 
import timeit

start = timeit.default_timer()
nA = 6.023e23
vol = 1e-15
X = [500,200,0,0]
c = [0,0,0]
c[0] = 2e-3
c[1] = 4e-2
c[2] = 0.1
V = [[-1,-1,1,0],[1,1,-1,0],[0,1,-1,1]]
a=[0,0,0]
t=0
tfinal = 30
tau = 0.1 

timeaxis = [0]
product1attime = [X[0]]
product2attime = [X[1]]
product3attime = [X[2]]
product4attime = [X[3]]

while t<tfinal:
    a[0] = c[0]*X[0]*X[1]
    a[1] = c[1]*X[2]
    a[2] = c[2]*X[2]
    p0 = np.random.poisson(a[0]*tau)
    p1 = np.random.poisson(a[1]*tau)
    p2 = np.random.poisson(a[2]*tau)
    p = [p0,p1,p2]
    for j in range(len(X)):
        X[j] = X[j] + p0*V[0][j] + p1 * V[1][j] + p2*V[2][j]
        
    t = t + tau    
    timeaxis.append(t)
    product1attime.append(X[0])
    product2attime.append(X[1])
    product3attime.append(X[2])
    product4attime.append(X[3])   
stop = timeit.default_timer()

print('Time: ', stop - start)    
plt.plot(timeaxis,product1attime, 'x')
plt.plot(timeaxis,product2attime, 'x')   
plt.plot(timeaxis,product3attime, 'x')   
plt.plot(timeaxis,product4attime, 'x')  
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.show()