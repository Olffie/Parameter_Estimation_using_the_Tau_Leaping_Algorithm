# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 11:48:37 2021

@author: Olaf
"""


"Tau leaping with linear parameter, chapter 3"

import numpy as np
from random import random
import matplotlib.pyplot as plt 
from Newton_Raphson_algorithm import findnewX, checkdifference



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
concentrationattime = [X.copy()]
product1attime = [X[0]]
product2attime = [X[1]]
product3attime = [X[2]]
product4attime = [X[3]]
C0 =  [X[0]]
C1 = [X[1]]
reactions0 = []

# the proces
while t<tfinal:
    c[0] = 2e-3 + 2e-4*t
    c[1] = 4e-2 + 4e-3 *t
    c[2] = 0.1 + 0.01*t
    a[0] = c[0]*X[0]*X[1]
    a[1] = c[1]*X[2]
    a[2] = c[2]*X[2]
    p0 = np.random.poisson(a[0]*tau)
    p1 = np.random.poisson(a[1]*tau)
    p2 = np.random.poisson(a[2]*tau)
    p = [p0,p1,p2]
    reactions0.append(p0)
    for j in range(len(X)):
        X[j] = X[j] + p0*V[0][j] + p1 * V[1][j] + p2*V[2][j]
        
    t = t + tau    
    timeaxis.append(t)
    product1attime.append(X[0])
    product2attime.append(X[1])
    product3attime.append(X[2])
    product4attime.append(X[3]) 
    C0.append(X[0]) 
    C1.append(X[1])
    concentrationattime.append(X.copy())
    
estimatorb0 =[]
estimatorc0 = []   
# find MLE
for i in range(2,len(C0)):
    
    xold = [10,10]
    xnew = [0.0001,0.0001]
    reaction0part = reactions0[0:i]
    while checkdifference(xold[0],xnew[0],xold[1],xnew[1],0.0000001) != True:
        xold = xnew
        xnew = findnewX(xold, reaction0part, C0, C1, tau)
    estimatorb0.append(xnew[0])
    estimatorc0.append(xnew[1])

upperboundc0 = []
lowerboundc0 = []

for i in range(2,len(estimatorc0)):
    
    sum1 = 0
    sum2 = 0
    sum3 = reactions0[0]/(estimatorb0[0]*0+estimatorc0[0])
    # make Fisher information matrix
    for j in range(1,i):

        sum1 = sum1 + reactions0[j]*(tau*j)**2 / (estimatorb0[j]*j*tau + estimatorc0[j])**2
        sum2 = sum2 + reactions0[j]*(tau*j) / (estimatorb0[j]*j*tau + estimatorc0[j])**2
        sum3 = sum3 + reactions0[j] / (estimatorb0[j]*j*tau + estimatorc0[j])**2
    I = [[sum1, sum2], [sum2, sum3]]
    factor = I[0][0]*I[1][1]  - I[1][0]*I[0][1] 
    I_inverse = [[I[1][1]/factor,-I[0][1]/factor], [-I[1][0]/factor,I[0][0]/factor]]
    lowerboundc0.append(estimatorc0[i-1] -  I_inverse[1][1]**0.5 * 1.96)
    upperboundc0.append(estimatorc0[i-1] +  I_inverse[1][1]**0.5 * 1.96)



plt.plot(timeaxis,product1attime, 'x')
plt.plot(timeaxis,product2attime, 'x')   
plt.plot(timeaxis,product3attime, 'x')   
plt.plot(timeaxis,product4attime, 'x')  
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.show()



print("a0 down", estimatorb0[298]-  I_inverse[0][0]**0.5 * 1.96)
print("a0 up", estimatorb0[298]+  I_inverse[0][0]**0.5 * 1.96)
print('a0 ', estimatorb0[298])



print("value b0:" , estimatorc0[296])
print("lowerbound:" , lowerboundc0[296])
print("upperbound: " , upperboundc0[296])
timeaxis.pop(300)
timeaxis.pop(0)
lowerboundc0.append(0)
lowerboundc0.append(0)
plt.plot(timeaxis,estimatorc0)
plt.plot(timeaxis,lowerboundc0)
plt.plot(timeaxis,[1e6 / (nA*vol)]*299)
plt.title("")
plt.show() 

