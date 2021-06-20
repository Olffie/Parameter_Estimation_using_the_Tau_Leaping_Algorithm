# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 10:52:40 2021

@author: Olaf

Determine the confidence interval
Chapter 5
"""
from sympy import symbols, Eq, solve
import numpy as np

def solve_aggregation(product,reactions,n,L,R,nonzeros):
    if (R-L)<0.0001:
        return (R+L)/2
    sum1 = calculate_sum1(reactions,n)
    sum2 = calculate_sum2(product,reactions,n,(L+R)/2,nonzeros)
    value = calculate_sum4(product,reactions,n,(L+R)/2,sum1,sum2)
    print(L, R, "value of middle point:" , value)
    if value<0: #het is een monotoon dalende functie
        return solve_aggregation(product,reactions,n,L,(L+R)/2,nonzeros)
    return solve_aggregation(product,reactions,n,(L+R)/2,R,nonzeros)


def calculate_sum1(reactions,n):
    sum_react = 0
    for i in range(n):#determine the total amount of reations
        for j in range(i,n):
            for k in range(len(reactions)):
                sum_react = sum_react + reactions[k][i][j]
    return sum_react

def calculate_sum2(product,reactions,n,alfa,nonzeros):
    sum_total= 0
    for k in range(len(reactions)): 
        for i in nonzeros[k]:
            for j in nonzeros[k]:
                if i<=j: #prevent doubles
                    if i==j:
                        sum_total =sum_total + ((i+1)*(1+j))**alfa * product[k][i]*product[k][j]/2          
                    else:  
                        sum_total =sum_total + ((i+1)*(1+j))**alfa * product[k][i]*product[k][j]           
    return sum_total

def calculate_sum3(product,reactions,n,alfa,i,j,sum1,sum2):
    sum_total = 0
    for k in range(len(reactions)):
        if ((product[k][i]!=0) & (product[k][j]!=0)):
            if i==j:
                sum_total = sum_total + reactions[k][i][j] - sum1/sum2 * ((i+1)*(1+j))**alfa * product[k][i]*product[k][j]/2
            else:
                sum_total = sum_total + reactions[k][i][j] - sum1/sum2 * ((i+1)*(1+j))**alfa * product[k][i]*product[k][j]
    return sum_total   

def calculate_sum4(product,reactions,n,alfa,sum1,sum2):
    sum_total= 0
    for i in range(n):
        for j in range(i,n):
            sum_total = sum_total + np.log((i+1)*(1+j))*calculate_sum3(product,reactions,n,alfa,i,j,sum1,sum2)
    return sum_total     

def determine_C(product,reactions,n,alfa,tau,nonzeros):
    sum1 = calculate_sum1(reactions,n)
    sum2 = calculate_sum2(product,reactions,n,alfa,nonzeros)
    return (sum1/sum2)

def Fisherinformationmatrix11(reactions,n,MLE_C):
    return calculate_sum1(reactions,n)/MLE_C**2

def Fisherinformationmatrix22(product,reactions,n,MLE_alfa, MLE_C):
    sum_total= 0
    for i in range(n): 
        for j in range(i,n):
            for k in range(len(reactions)): 
                if ((product[k][i]!=0) & (product[k][j]!=0)):#bc the last one does not reaciton anymore
                    if i==j:
                        sum_total =sum_total + MLE_C* (np.log((i+1)*(1+j)))**2 * ((i+1)*(1+j))**MLE_alfa * product[k][i]*product[k][j]/2           
                    else:
                        sum_total =sum_total + MLE_C* (np.log((i+1)*(1+j)))**2 * ((i+1)*(1+j))**MLE_alfa * product[k][i]*product[k][j]    
    return sum_total

def Fisherinformationmatrix12(product,reactions,n,MLE_alfa):
    sum_total= 0
    for i in range(n): 
        for j in range(i,n):
            for k in range(len(reactions)):
                if ((product[k][i]!=0) & (product[k][j]!=0)):
                    if i==j: 
                        sum_total =sum_total + np.log((i+1)*(1+j)) * ((i+1)*(1+j))**MLE_alfa * product[k][i]*product[k][j]/2
                    else:
                        sum_total =sum_total + np.log((i+1)*(1+j)) * ((i+1)*(1+j))**MLE_alfa * product[k][i]*product[k][j]
    return sum_total

def InverseFishermatrix(product,reactions,n,MLE_alfa, MLE_C):
    sum1 = Fisherinformationmatrix11(reactions,n,MLE_C)
    sum2 = Fisherinformationmatrix12(product,reactions,n,MLE_alfa)
    sum3 = Fisherinformationmatrix22(product,reactions,n,MLE_alfa, MLE_C)
    I = [[sum1, sum2], [sum2, sum3]]
    factor = I[0][0]*I[1][1]  - I[1][0]*I[0][1] 
    I_inverse = [[I[1][1]/factor,-I[0][1]/factor], [-I[1][0]/factor,I[0][0]/factor]]
    return I_inverse

def Fisherinformation_C(product,reactions, n, MLE_C, MLE_alfa):
    return InverseFishermatrix(product,reactions,n,MLE_alfa, MLE_C)[0][0]

def Fisherinformation_alfa(product,reactions, n, MLE_C, MLE_alfa):
    return InverseFishermatrix(product,reactions,n,MLE_alfa, MLE_C)[1][1]