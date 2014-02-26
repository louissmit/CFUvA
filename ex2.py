# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/henning/.spyder2/.temp.py
"""

import numpy
import random
import matplotlib.pyplot as plt


def computeMCValuation(k=99.0, s0=100.0, r=0.06, v=0.2, T=1, M=100000):
    sols = numpy.zeros((M,1))
    
    for i in xrange(0,M):
        Z = random.gauss(0,1)
        sn = s0*numpy.exp((r-0.5*pow(v,2))*T + v*numpy.sqrt(T)*Z)
        
        sols[i] = max(0,sn-k)*numpy.exp(-r*T)
        
    
    return numpy.mean(sols), numpy.sqrt(numpy.var(sols))/numpy.sqrt(M)
    
    
def computeDelta(payoff = 'EuropeanC',k=99.0, s0=100.0, r=0.06, v=0.2, T=1, M=100000, bump = 0.1, sameSeed = True):
    
    shocks = numpy.zeros((M,1))
    
    for i in xrange(0,M):
        Z = random.gauss(0,1)
        shocks[i] = numpy.exp((r-0.5*pow(v,2))*T + v*numpy.sqrt(T)*Z)
        
    s = numpy.array(shocks)*s0

    if not sameSeed:
        for i in xrange(0,M):
            Z = random.gauss(0,1)
            shocks[i] = numpy.exp((r-0.5*pow(v,2))*T + v*numpy.sqrt(T)*Z)
    
    
    s_bumped = numpy.array(shocks)*(s0+bump)
    
    #computing payoffs
    if payoff == 'EuropeanC':
        s_payoff = [max(x-k,0) for x in s]
        sb_payoff = [max(x-k,0) for x in s_bumped]
    elif payoff == 'Digital':
        s_payoff = [x>k for x in s]
        sb_payoff = [x>k for x in s_bumped]
        
    
    return (numpy.mean(sb_payoff)-numpy.mean(s_payoff))/bump
    
    
def computeDeltaApprox(k=99.0, s0=100.0, r=0.06, v=0.2, T=1, M=100000):
    s = numpy.zeros((M,1))    
    
    for i in xrange(0,M):
        Z = random.gauss(0,1)
        s[i] = s0*numpy.exp((r-0.5*pow(v,2))*T + v*numpy.sqrt(T)*Z)
    
    
    sigmoid = lambda x,k: 1/(1+numpy.exp(-(x-k)))
    sigmoid_d = lambda x,k: sigmoid(x,k)*(1-sigmoid(x,k))
    
    sigmoid_d2 = lambda x,k: numpy.abs(x-k)<0.5
    
    deltas = [sigmoid_d2(x,k)*x/s0 for x in s]    
    
    return numpy.mean(deltas)
        

def computeDeltaLikelihood(k=99.0, s0=100.0, r=0.06, v=0.2, T=1, M=100000):
    shocks = numpy.zeros((M,1))
    shocks = [random.gauss(0,1) for x in shocks]
    
    
    s_t = lambda Z: s0*numpy.exp((r-0.5*pow(v,2))*T + v*numpy.sqrt(T)*Z)
    weird_function = lambda Z: numpy.exp(-r*T)*(s_t(Z)>k)*Z/(v*s0*numpy.sqrt(T))
    
    deltas = [weird_function(Z) for Z in shocks]
    
    return numpy.mean(deltas)
    
    
def plotValuationConvergence(n=10, ms = numpy.linspace(50,100000,15)):
    X=[]
    Y=[]
    
    for m in ms:
        for j in xrange(0,n):
            #print int(m)
            X.append(int(m))
            Y.append(computeMCValuation(M=int(m))[0])
            
    print X,Y
    plt.scatter(X,Y)
    
def varianceAntithetic(k=99.0, s0=100.0, r=0.06, v=0.2, T=1, M=100, meanPayoff=True):
    shocks = numpy.zeros((M,1))
    shocks = [random.gauss(0,1) for x in shocks]
    #shocks = [max(x,-x) for x in shocks]
    shocks2= numpy.array(shocks)*(-1)
    
    #shocks = numpy.concatenate((shocks,shocks2),axis=1)
    
    
    s_t = lambda Z: s0*numpy.exp((r-0.5*pow(v,2))*T + v*numpy.sqrt(T)*Z)
    
    sn = [s_t(Z) for Z in shocks]
    sn2 = [s_t(Z) for Z in shocks2]
    
    if not meanPayoff:
        sn = (numpy.array(sn)+numpy.array(sn2))/2
    
    if meanPayoff:
        payoffs = [max(0,x-k)*numpy.exp(-r*T) for x in sn]
        payoffs2 = [max(0,x-k)*numpy.exp(-r*T) for x in sn2]
    
        payoffs = (numpy.array(payoffs)+numpy.array(payoffs2))/2
    else:
        payoffs = [max(0,x-k)*numpy.exp(-r*T) for x in sn]
    
    #return sn
    return numpy.mean(payoffs), numpy.sqrt(numpy.var(payoffs))/numpy.sqrt(M)
    
    
def plotAntitheticShit():
    total_M = 10000
    m = numpy.linspace(50,5000,50)
    
    X = []
    Y2 = []
    Y = []
    Y3 = []
    
    for mi in m:
        results1 = []
        results2 = []
        results3 = []
        for i in xrange(0,int(total_M/mi)):
            results1.append(varianceAntithetic(M=int(mi/2))[1])
            results3.append(varianceAntithetic(M=int(mi/2), meanPayoff = False)[1])
            results2.append(computeMCValuation(M=int(mi))[1])
        Y.append(numpy.mean(results1))
        Y2.append(numpy.mean(results2))
        Y3.append(numpy.mean(results3))
        X.append(mi)
        
    plt.plot(X,Y,label='Antithetic Payoffmean')
    plt.plot(X,Y,label='Antithetic Stockmean')
    plt.plot(X,Y2,label='Vanilla')
    plt.legend()
    plt.show()

print computeMCValuation()