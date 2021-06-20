# -*- coding: utf-8 -*-
"""
Created on Sat May 15 13:23:26 2021

@author: Olaf

Plot of the data, chapter 6
"""
import numpy as np
from Aggregation_solve_MLE_realdata import solve_aggregation, determine_C, Fisherinformation_alfa, Fisherinformation_C
from scipy.stats import poisson
import matplotlib.pyplot as plt 

def getconcentration(product,i):
    concentrationx = []
    for j in range(len(product)):
        concentrationx.append(product[j][i])
    return concentrationx


def findreaction(low,high,reactions,time):
    reactionsattime = reactions[time]
    for i in range(len(reactionsattime)):
        reaction = reactionsattime[i]
        if (reaction[0] == high) and (reaction[1] == low):
            return reaction[2]
    return 0

def highestindex(X):
    i = len(X)-1
    while i>0:
        if X[i]!=0:
            return i+1
        i-=1
timeaxis = [0]

with open('sizes.txt', 'r') as sizes:
    lines = sizes.readlines()    
X = [[0]*2000] #add starting situation
X[0][0] = 2000
for i in range(len(lines)):
    con = lines[i].split()
    lengthcon = len(con)
    for i in range(lengthcon,2000):
        con.append(0)
    #if(len(X))<89:    
    X.append( con)

for i in range(len(X)): 
    for j in range(len(X[i])):
        X[i][j] = int(X[i][j])
X.pop(0)        
nonzeros = []
existing = set()
for k in range(len(X)):
    nonzero = []
    for i in range(len(X[k])):
        if X[k][i]!=0:
            nonzero.append(i)
            existing.add(i)
    nonzeros.append((nonzero))     
reactions = []
for k in range(1,2000):
    try:
        filename = str(k) + ".txt"
        R = []
        with open(filename, 'r') as reaction:
            line = reaction.readlines()
            for l in range(len(line)):
                r = list(map(int, line[l].split()))
                r[0] = r[0]-1 #fixing the indices
                r[1] = r[1]-1
                R.append(r)
        reactions.append(R)
        timeaxis.append(k)
        
    except:
         tau = 1

tau = 1
maximumsize = []
for i in range(len(X)):
    maximumsize.append(highestindex(X[i]))

print(len(timeaxis))
#plt.plot(timeaxis,getconcentration(X,0), 'x')    
#plt.plot(timeaxis,getconcentration(X,1), 'x')
#plt.plot(timeaxis,getconcentration(X,3), 'x')  
#plt.plot(timeaxis,getconcentration(X,5), 'x')
#plt.plot(timeaxis,maximumsize,'x')
plt.yscale("log")
plt.xlabel("Time scaled")
plt.ylabel("Size largest cluster")
plt.show()


len(X)

MLE_alfa = solve_aggregation(X,reactions,2000,0,2,nonzeros,existing)
MLE_C = determine_C(X,reactions,2000,MLE_alfa,tau,nonzeros,existing) 
print(MLE_alfa)
print(MLE_C)  


FisherC = Fisherinformation_C(X,reactions, 2000, MLE_C, MLE_alfa,existing)
R_side_C = MLE_C/tau + 1.96*FisherC**0.5/tau
L_side_C = MLE_C/tau - 1.96*FisherC**0.5/tau
Fisheralfa = Fisherinformation_alfa(X,reactions, 2000, MLE_C, MLE_alfa,existing)

print("Confidence interval C:" , L_side_C, R_side_C)
R_side_alfa = MLE_alfa + 1.96*Fisheralfa**0.5 ##I have to use to 2D version
L_side_alfa = MLE_alfa - 1.96*Fisheralfa**0.5
print("Confidence interval alfa:" , L_side_alfa, R_side_alfa)


  