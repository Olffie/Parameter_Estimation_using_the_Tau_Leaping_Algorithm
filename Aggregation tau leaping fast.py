# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 16:41:14 2021

@author: Olaf
Chapter 5
"Aggregation Tau leaping (fast implementation)"
For scaled time and non-scaled time
"""


import numpy as np
from random import random
import matplotlib.pyplot as plt 
from Aggregation_solve_MLE import solve_aggregation, determine_C, Fisherinformation_alfa, Fisherinformation_C

def Cloning(li1):
    li_copy = []
    for i in range(len(li1)):
        li_copy.append(li1[i].copy())
    return li_copy
def getconcentration(product,i):
    concentrationx = []
    for j in range(len(product)):
        concentrationx.append(product[j][i])
    return concentrationx

def summie(X):
    som = 0
    for i in range(len(X)):
        som = som + (i+1)*X[i]
    return som
def highestindex(X):
    i = len(X)-1
    while i>0:
        if X[i]!=0:
            return i+1
        i-=1
        
erin_alfa = 0
erin_C = 0;


n = 2000
X = []
C = []
c = []
maximumsize = []
alfa = 0.481597900390625
constant = 2.4732524807359215e-05
for i in range(n): #initial conditions
    X = [0]*n
    X[0] = 2000
    #X.append(np.math.ceil(1000/(i+1)))
for i in range(n):
    c.append(constant) #c/2 case in poisson     
for i in range(n):
    C.append(c)
        


t=0
tfinal = 18000
tau = 0.5 #0.02 van maken

timeaxis = [0]
product = [X.copy()]
reactions = []


    
while t<tfinal:
    a = np.zeros((n,n), dtype = 'f2')
    p = np.zeros((n,n), dtype = 'f2')
    maximumsize.append(highestindex(X))
    nonzero = []
    
    for i in range(n):
        if X[i]!= 0:
            nonzero.append(i)
    for i in nonzero:
        for j in nonzero:
            if i<=j:
                
                if i==j:
                    a[i][j] = 0.5*C[i][j]*X[i]*X[j]*((i+1)*(j+1))**alfa 
                else: 
                    a[i][j] = C[i][j]*X[i]*X[j]*((i+1)*(j+1))**alfa 
                 
    for i in nonzero:
        for j in nonzero:
            if a[i][j]!= 0:
                p[i][j] = np.random.poisson(a[i][j]*tau)               
    for i in nonzero:
        for j in nonzero:
            if p[i][j]!= 0 & i<=j:
                if i!=j:
                    if (X[i] - p[i][j] >= 0) & (X[j] - p[i][j] >= 0):
                        X[i] = X[i] - p[i][j]
                        X[j] = X[j] - p[i][j]
        
                        if i+j+1<n:
                            X[i+j+1] = X[i+j+1] + p[i][j]
                else:
                    if (X[i] - 2*p[i][j] >= 0):
                        X[i] = X[i] - p[i][j]
                        X[j] = X[j] - p[i][j]
            
                        if i+j+1<n:
                            X[i+j+1] = X[i+j+1] + p[i][j]

    t = t + tau 
    # in case of scaled time
    timeaxis.append(n-sum(X))
    # in case of lineair time
    # timeaxis.append(t)
    print("time: " , t)
    product.append(X.copy())
        
    #reactions.append(Cloning(p))


maximumsize.append(highestindex(X))
plt.plot(timeaxis,getconcentration(product,0), 'x')    
plt.plot(timeaxis,getconcentration(product,1), 'x')
plt.plot(timeaxis,getconcentration(product,3), 'x')  
plt.plot(timeaxis,getconcentration(product,5), 'x')  
#plt.plot(timeaxis,maximumsize, 'x')

plt.yscale("log")
plt.xlabel("Time scaled")
plt.ylabel("Size largest cluster")
plt.show()








nonzeros = []
for k in range(len(product)):
    nonzero = []
    for i in range(len(product[k])):
        if product[k][i]!=0:
            nonzero.append(i)
    nonzeros.append((nonzero))




MLE_alfa = solve_aggregation(product,reactions,n,0,2,nonzeros)
MLE_C = determine_C(product,reactions,n,MLE_alfa,tau,nonzeros) ##tau is gone
FisherC = Fisherinformation_C(product,reactions, n, MLE_C, MLE_alfa) 
Fisheralfa = Fisherinformation_alfa(product,reactions, n, MLE_C, MLE_alfa)
print("MLE_C:" , MLE_C/tau)
print("MLE_alfa:" , MLE_alfa)
R_side_C = MLE_C/tau + 1.96*FisherC**0.5/tau
L_side_C = MLE_C/tau - 1.96*FisherC**0.5/tau
R_side_alfa = MLE_alfa + 1.96*Fisheralfa**0.5 ##I have to use to 2D version
L_side_alfa = MLE_alfa - 1.96*Fisheralfa**0.5
print("Confidence interval alfa:" , L_side_alfa, R_side_alfa)
print("Confidence interval C:" , L_side_C, R_side_C) 
if L_side_alfa < alfa <R_side_alfa:
    erin_alfa += 1
if L_side_C < constant<R_side_C:
    erin_C += 1

     


