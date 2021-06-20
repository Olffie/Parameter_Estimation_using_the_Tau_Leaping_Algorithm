# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 11:53:06 2021

@author: Olaf


"""
"Newton-Raphson algorithm for the lineair rate constant"

def checkdifference1(b0,b1,a):
    if(abs(b0-b1)>a):
        return False
    return True

def checkdifference(b0,b1,c0,c1,a):
    if(checkdifference1(b0,b1,a) & checkdifference1(c0,c1,a)):
        return True
    return False

def jacobian(X,b0,c0,tau):
    n = len(X)
    a11 = 0
    a12 = 0
    a21 = 0
    a22 = 0
    for i in range(n):
        a11 = a11 + -X[i]*(tau*i) * (b0*(tau*i)+c0)**-2
        a12 = a12 + -X[i] * (b0*(tau*i)+c0)**-2
        a21 = a21 + -X[i] *(tau*i)**2 * (b0*(tau*i)+c0)**-2
        a22 = a22 + -X[i] *(tau*i) * (b0*(tau*i)+c0)**-2
    return [[a21,a22],[a11,a12]] #rijen verkeerd opgeschreven

def inverse2x2Jacobian(J):
    factor = J[0][0]*J[1][1]  - J[1][0]*J[0][1] 
    return [[J[1][1]/factor,-J[0][1]/factor], [-J[1][0]/factor,J[0][0]/factor]]

def findnewX(x,X,C0,C1,tau):

    xnew = [0,0]
    fx = [0,0]
    b0 = x[0]
    c0 = x[1]
    J = jacobian(X,b0,c0,tau)

    Jinverse = inverse2x2Jacobian(J)

    for i in range(len(X)):

        fx[0] = fx[0] + X[i]*(tau*i) / (b0*(tau*i)+c0) -C0[i]*C1[i]*tau*(tau*i)

        fx[1] = fx[1] + X[i] / (b0*(tau*i)+c0) -C0[i]*C1[i]*tau
        

    for i in range(2):
        xnew[i] = x[i] - (Jinverse[i][0] * fx[0] + Jinverse[i][1]*fx[1])
    return xnew
        
        
    