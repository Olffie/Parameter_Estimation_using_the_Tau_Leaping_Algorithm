# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 12:36:51 2021

@author: Olaf
Chapter 3/4
"""

"Tau leaping with (Estimation one reaction) -> PD/normal"

import numpy as np
from random import random
import matplotlib.pyplot as plt 
from PD_distribution import findconstants, FishermatrixinversePD

nA = 6.023e23
vol = 1e-15
c = [0,0,0]
c[0] = 2e-3
c[1] = 4e-2
c[2] = 0.1
V = [[-1,-1,1,0],[1,1,-1,0],[0,1,-1,1]]


tau = 0.1
aantal = 100
erin = 0

X = [500,200,0,0]
a=[0,0,0]
t=0
tfinal = 30
timeaxis = [0]
concentrationattime = [X.copy()]
productattime = [X[2]]
substrateattime =  [X[0]]
reactions0 = []

X_var =[]
Y_var = []
Z_var = []



while t<tfinal:
    X_var.append(X[0]*X[1]*tau)
    Y_var.append(X[2]*tau)
    a[0] = c[0]*X[0]*X[1]
    a[1] = c[1]*X[2]
    a[2] = c[2]*X[2]
    p0 = np.random.poisson(a[0]*tau)
    p1 = np.random.poisson(a[1]*tau)
    p2 = np.random.poisson(a[2]*tau)
    Z_var.append(p0-p1)
    p = [p0,p1,p2]
    reactions0.append(p0)
    for j in range(len(X)):
        X[j] = X[j] + p0*V[0][j] + p1 * V[1][j] + p2*V[2][j]
        
    t = t + tau    
    timeaxis.append(t)
    productattime.append(X[2])
    substrateattime.append(X[0])  
    concentrationattime.append(X.copy())
    

  
Con = findconstants(X_var, Y_var, Z_var)
print(Con)
I_inverse = FishermatrixinversePD(X_var,Y_var,Z_var,Con[0],Con[1])
print(Con[0]-  I_inverse[0][0]**0.5 * 1.96, Con[0]+  I_inverse[0][0]**0.5 * 1.96 )
print(Con[1]-  I_inverse[1][1]**0.5 * 1.96, Con[1]+  I_inverse[1][1]**0.5 * 1.96 )  

'''
plt.plot(timeaxis,productattime, 'x')    
plt.plot(timeaxis,substrateattime, 'x')   
plt.xlim(0, 60)
plt.ylim(0,1000)
plt.title("Tau leaping")
#print(erin)
'''
#PD difference
Con = findconstants(X_var, Y_var, Z_var)
I_inverse = FishermatrixinversePD(X_var,Y_var,Z_var,Con[0],Con[1])
print(Con[0]-  I_inverse[0][0]**0.5 * 1.96, Con[0]+  I_inverse[0][0]**0.5 * 1.96 )
print(Con[1]-  I_inverse[1][1]**0.5 * 1.96, Con[1]+  I_inverse[1][1]**0.5 * 1.96 )

print(erin)

#normal
sum1 = reactions0[0]
sum2 = concentrationattime[0][0] * concentrationattime[0][1]*tau

for i in range(1,len(productattime)-1):
    sum1 = sum1 + reactions0[i]
    sum2 = sum2 + concentrationattime[i][0] * concentrationattime[i][1]*tau

estimatorc0 = (sum1/sum2)
variance2 = (1.96*(-(-estimatorc0**-2 * sum1))**-0.5 ) #toevoeging 0.001
upperbound = (estimatorc0 + variance2)
lowerbound = (estimatorc0 - variance2)
print(estimatorc0 -variance2, estimatorc0 + variance2)


'''
timeaxis.pop(300)
timeaxis.pop(0)
plt.plot(timeaxis,estimatorc0)
plt.plot(timeaxis, upperbound)
plt.plot(timeaxis,lowerbound)
plt.plot(timeaxis,[0.002]*299)
plt.title("Estimation parameter with 95% confidence interval")
plt.show()
'''

