# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 15:22:01 2021

@author: Olaf

Chapter 2
"""
import numpy as np
from random import random
import matplotlib.pyplot as plt 
import timeit

def findindex(random, prop):
    for i in range(len(prop)):
        if random<prop[i]:
            return i
    return -1
  
def determine_index_v(xprev,xnew, V):
    for i in range(len(V)):
        if np.add(xprev,V[i]) == xnew:
            return i
    return -1

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
timeaxis = [0]
product1attime = [X[0]]
product2attime = [X[1]]
product3attime = [X[2]]
product4attime = [X[3]]

concentrationprojectory = [X]
while t<tfinal:
    a[0] = c[0]*X[0]*X[1]
    a[1] = c[1]*X[2]
    a[2] = c[2]*X[2]
    asum = sum(a)
    prop = np.cumsum(a)/asum
    
    rand = random()
    index = findindex(rand, prop)
    tau = np.log(1/rand)/asum
    X = np.add(X,V[index])
    t = t + tau
    timeaxis.append(t)
    product1attime.append(X[0])
    product2attime.append(X[1])
    product3attime.append(X[2])
    product4attime.append(X[3])

    concentrationprojectory.append(X)



#plot the graph
plt.plot(timeaxis,product1attime, 'x', label='S_1')
plt.plot(timeaxis,product2attime, 'x')   
plt.plot(timeaxis,product3attime, 'x')   
plt.plot(timeaxis,product4attime, 'x')       
plt.xlabel("Time")
plt.ylabel("Concentration")

#plt.xlim(12, 17)
#plt.ylim(80,130)
plt.show()

    
    
    
