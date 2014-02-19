# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/henning/.spyder2/.temp.py
"""

import numpy
import random


def computeMCValuation(k=99, s0=100, r=0.06, v=0.2, T=1, M=100000):
    sols = numpy.zeros((M,1))
    
    for i in xrange(0,M):
        Z = random.gauss(0,1)
        sn = s0*numpy.exp((r-0.5*pow(v,2))*T + v*sqrt(T)*Z)
        
        sols[i] = max(0,sn-k)*numpy.exp(-r*T)
        
    
    return numpy.mean(sols), numpy.var(sols)/sqrt(M)