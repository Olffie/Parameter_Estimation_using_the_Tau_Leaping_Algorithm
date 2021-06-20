# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 11:07:59 2021

@author: Olaf
Chapter 4
"""


from scipy.special import iv
from sympy import solve, Eq
from sympy.abc import x,y

from numpy import math

def findvalue(X,Y,Z,theta2):
    sumX = sum(X)
    sumY = sum(Y)
    sumZ = sum(Z)
    theta1 = (sumZ + theta2*sumY)/sumX
    besselsum  = 0
    for i in range(len(Z)):
        if Y[i]==0:
            Y[i] = 1e-16
        sqrtex = (X[i]*Y[i]*theta1*theta2)**0.5
        besselsum = besselsum + X[i]*Y[i]*theta2/sqrtex * iv(Z[i]+1,2*sqrtex)/iv(Z[i],2*sqrtex)
    value = -sumX + sumZ/theta1 + besselsum
    return value

def findmax(X,Y,Z,theta2):
    sumX = sum(X)
    sumY = sum(Y)
    sumZ = sum(Z)
    theta1 = (sumZ + theta2*sumY)/sumX
    
    sum1  = 0
    for i in range(len(Z)):
        sqrtex = (X[i]*Y[i]*theta1*theta2)**0.5
        sum1 = sum1 + -theta1*X[i]-theta2*Y[i]+Z[i]/2 *math.log(X[i]*theta1/(Y[i]*theta2)) + math.log(iv(Z[i],2*sqrtex))
    return sum1

def solve_eq(X,Y,Z,step, start):
    if step<0.000000001:
        return start
    best = 1e20
    best_start  = 0
    for i in range(-10,10):
        
        value = findvalue(X,Y,Z,start + i*step)
        if abs(value-0)<best:
            best = abs(value-0)
            best_start = start + i*step        
    return solve_eq(X,Y,Z,step/10,best_start)
            
def findconstants(X,Y,Z):
    theta2 = solve_eq(X,Y,Z,0.1,1.1111111)
    sumX = sum(X)
    sumY = sum(Y)
    sumZ = sum(Z)
    theta1 = (sumZ + theta2*sumY)/sumX
    return [theta1, theta2]
    
def Fishermatrix(X,Y,Z,c1,c2):
    sum11 = 0
    sum12 = 0
    sum22 = 0
    for i in range(1,len(Z)):
        xi = 2*(X[i]*Y[i]*c1*c2)**0.5
        bes0 = iv(Z[i],xi)
        bes1 = iv(Z[i]+1,xi)
        bes2 = iv(Z[i]+2,xi)
        sum11 = sum11 + Z[i]/c1 + 0.5*c1**-1.5 * (X[i]*Y[i]*c2)**0.5 * bes1/bes0 + X[i]*Y[i]*c2/c1*(xi**-1 * bes1*bes0 + bes2*bes0 - bes0**2)/bes0**2 
        sum22 = sum22  + 0.5*c2**-1.5 * (X[i]*Y[i]*c1)**0.5 * bes1/bes0 + X[i]*Y[i]*c1/c2*(xi**-1 * bes1*bes0 + bes2*bes0 - bes0**2)/bes0**2 
        sum12 = sum12 - 0.5*xi * bes1/bes0 + X[i]*Y[i]*(xi**-1 * bes1*bes0 + bes2*bes0 - bes0**2)/bes0**2 
    I = [[sum11, sum12], [sum12, sum22]]
    return I

def FishermatrixinversePD(X,Y,Z,c1,c2):
    I = Fishermatrix(X,Y,Z,c1,c2)
    factor = I[0][0]*I[1][1]  - I[1][0]*I[0][1] 
    I_inverse = [[I[1][1]/factor,-I[0][1]/factor], [-I[1][0]/factor,I[0][0]/factor]]
    return I_inverse
    
