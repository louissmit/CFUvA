# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 11:25:55 2014

@author: henning
"""

import numpy

def LETSDOIT(I=100,N=100, r=0.06, v=0.2, s0 = 100.0, k = 99.0, T=1.0):
    big_mat = numpy.zeros((I+1,N+1))
    int_incrs = 2*s0/I
    
    delta_t = T/N    
    
    for i in xrange(0,I+1):
        big_mat[i,0]=max(int_incrs*i-k,0)
        
    
    a1=0
    a0=1
    a_1=0    
    
    A = numpy.zeros((N+1,N+1))
    #FTCS:
    for i in xrange(0,N+1):
        A[i,i]=1
    Abar = numpy.linalg.inv(A)
    
    alpha = (r-pow(v,2)/2)*(delta_t/(2*int_incrs))
    beta = (pow(v,2)/2)*(delta_t/pow(int_incrs,2))
    gamma = r*delta_t
    
    b1 = (alpha + beta)
    b0 = (1-2*beta-gamma)
    b_1 = (beta - alpha)
    
    print b1,b0,b_1
    
    for i in xrange(1,N+1):
        #calculate c
        c = numpy.zeros((I+1,1))
        #not sure if c[0] and c[I] make sense
        c[0] = b1*big_mat[1,i-1]+b0*big_mat[0,i-1] #+b_1*big_mat[i-1,0-1]-a_1*big_mat[i,0-1]
        #c[I] = b1*big_mat[i-1,j+1]+b0*big_mat[i-1,j]+b_1*big_mat[i-1,j-1]-a1*big_mat[i,j+1]
        c[I] = b0*big_mat[N,i-1]+b_1*big_mat[N-1,i-1]
        for j in xrange(1,I):
            c[j] = b1*big_mat[j+1,i-1]+b0*big_mat[j,i-1]+b_1*big_mat[j-1,i-1]
        big_mat[:,i] = numpy.dot(Abar,c).T
        
    
    return big_mat
    