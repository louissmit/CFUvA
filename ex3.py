# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 11:25:55 2014

@author: henning
"""
import matplotlib.pyplot as plt
import numpy
import pprint
from blackscholes import BS

def FD(I=100,N=100, r=0.06, v=0.2, s0 = 100.0, k = 99.0, T=1.0, type='ftcs'):
    V = numpy.zeros((I+1,N+1))
    M1 = -2.0
    M2 = 2.0
    
    l = numpy.log(s0)
    
    print 'smin', numpy.exp(l+M1), ' smax',numpy.exp(l+M2)

    V[:,0] = numpy.linspace(M1+l, l+M2, I+1)

    s0s = numpy.exp(numpy.linspace(M1+l, l+M2, I+1))    
    
    delta_x = V[2,0] - V[1,0]

    print delta_x    
    
    for i in xrange(0,I+1):
        V[i,0]=max(numpy.exp(V[i,0])-k,0)

    
    delta_t = T/N


    alpha = (r-pow(v,2)/2)*(delta_t/(2*delta_x))
    beta = (pow(v,2)/2)*(delta_t/pow(delta_x,2))
    gamma = r*delta_t

    if type == 'ftcs':
        a1=0
        a0=1
        a_1=0
    else:
        a1 = -alpha - beta/2
        a0 = beta + gamma/2 + 1
        a_1 = alpha - beta/2

    A = numpy.zeros((I+1, I+1))
    #FTCS:
    for i in xrange(0,I+1):
        A[i,i]=a0
        if i < I:
            A[i+1,i]=a1
            A[i,i+1]=a_1
            
    Abar = numpy.linalg.inv(A)
    

    

    if type == 'ftcs':
        b1 = (alpha + beta)
        b0 = (1-2*beta-gamma)
        b_1 = (beta - alpha)
    else:
        b1 = -a1
        b0 = -beta - gamma/2 + 1
        b_1 = -a_1

    print 'a', a1, a0, a_1
    print 'b', b1, b0, b_1

    for n in xrange(0,N):
        #calculate c
        c = numpy.zeros((I+1,1))
        #not sure if c[0] and c[I] make sense
        if type=='ftcs':
            c[0] = b1*V[2,n] + b0*V[1,n] + b_1*V[0, n] - a_1*V[0,n+1]
            c[I] = b0*V[I,n]+ b_1*V[I-1,n] + b1*(V[I,n]+delta_x*numpy.exp(l+M2))#- a1*V[I+1, n+1] + b1*V[I+1,n]
        else:
            c[0] = b1*V[2,n] + b0*V[1,n] + b_1*V[0, n] - a_1*V[0,n+1]
            c[I] = b0*V[I,n]+ beta*V[I-1,n] + b1*(V[I,n]+delta_x*numpy.exp(l+M2))

        for j in xrange(1,I):
            c[j] = b1*V[j+1,n] + b0*V[j,n] + b_1*V[j-1,n]

        V[:,n+1] = numpy.dot(Abar,c).T

    return V, s0s
    
    
def getError(N1=1000,I1=1000):
    mat,s = FD(N=N1, I=I1, type='cn')
    x = [BS(s0=s0)[0] for s0 in s]
    
    s = s[0:-1]
    x = x[0:-1]
    
    y = mat[0:-1,-1]
    
    print len(x)
    print len(y)
    print len(s)
    
    plt.plot(s,x,color='blue')
    plt.plot(s,y,color='red')
    

M,s = FD(type='cn')
getError()
