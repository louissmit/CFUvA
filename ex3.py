# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 11:25:55 2014

@author: henning
"""
import matplotlib.pyplot as plt
import numpy
import pprint

def LETSDOIT(I=100,N=50, r=0.06, v=0.2, s0 = 100.0, k = 99.0, T=1.0):
    V = numpy.zeros((I+2,N+1))
    M1 = -10.0
    M2 = 10.0
    
    V[1:,0] = numpy.linspace(M1, M2, I+1)
 
    int_incrs = (M2 - M1) / I
     
    delta_x = numpy.exp(0.2)
    delta_t = -T/N

    #for i in range_x:
    #    V[i,0]=i

    for n in xrange(0,N):
        V[I+1,n] = 1


    a1=0
    a0=1
    a_1=0

    A = numpy.zeros((I, I))
    #FTCS:
    for i in xrange(0,I):
        A[i,i]=1
    Abar = numpy.linalg.inv(A)

    alpha = (r-pow(v,2)/2)*(delta_t/(2*delta_x))
    beta = (pow(v,2)/2)*(delta_t/pow(delta_x,2))
    gamma = r*delta_t

    b1 = (alpha + beta)
    b0 = (1-2*beta-gamma)
    b_1 = (beta - alpha)
    
    print b1, b0, b_1


    for n in xrange(0,N):
        #calculate c
        c = numpy.zeros((I,1))
        #not sure if c[0] and c[I] make sense
        c[0] = b1*V[2,n]+b0*V[1,n] + b_1*V[0, n] #+b_1*big_mat[i-1,0-1]-a_1*big_mat[i,0-1]
        #c[I] = b1*big_mat[i-1,j+1]+b0*big_mat[i-1,j]+b_1*big_mat[i-1,j-1]-a1*big_mat[i,j+1]
        c[I-1] = b1*V[I+1,n]+b0*V[I,n]+b_1*V[I-1,n] #- a1*big_mat[I, n+1]
        for j in xrange(1,I-1):
            
            c[j] = b1*V[j+1,n]+b0*V[j,n]+b_1*V[j-1,n]
            
        V[0:I ,n+1] = numpy.dot(Abar,c).T
        #return c, numpy.dot(Abar,c)
        

    # plt.pcolormesh(big_mat)
    # plt.show()

    return V

M = LETSDOIT()
