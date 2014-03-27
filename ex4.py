# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 14:50:42 2014

@author: henning
"""

import numpy as np

s0 = 100.0
k = 80.0
T = 1.0
sig = 0.2
r = 0.06
a = np.log(s0/k) + r*T - 12*np.sqrt(pow(sig,2)*T)
b = np.log(s0/k) + r*T + 12*np.sqrt(pow(sig,2)*T)

X = np.log(s0/k)

#not sure if t = T
phi = lambda u: np.exp(1j*u*(r-pow(sig,2)/2)*T - pow(sig,2)*T*pow(u,2)/2)



def F(n):
    return 2/(b-a)*np.real(phi(n*np.pi/(b-a))*np.exp(-1j*n*a*X/(b-a)))
    
def G(n):
    weirdX = 1/(1+pow((n*np.pi)/(b-a),2))*(np.cos((n*np.pi)*np.exp(b)) - 
        np.cos((n*np.pi*a/(b-a))) + n*np.pi/(b-a)*np.sin(n*np.pi)*np.exp(b)
        -n*np.pi/(b-a)*np.sin(n*np.pi*a/(b-a)))
        
    #weirdV = 0
    if n == 0:
        weirdV = b
    else:
        weirdV = (b-a)/(n*np.pi)*(np.sin(k*np.pi) - np.sin(k*np.pi*a/(b-a)))
        
    return k*(weirdX-weirdV)*2/(b-a)
    
def cos_shit(n=124):
    summ = 0
    for i in xrange(0,n):
        summ += G(i)*F(i)
    return np.exp(-r*T)*(b-a)/2*summ