# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/henning/.spyder2/.temp.py
"""

import numpy
import scipy
import scipy.stats
import random
import matplotlib.pyplot as plt


def computeMCValuation(k=99.0, s0=100.0, r=0.06, v=0.2, T=1, M=100000):
    sols = numpy.zeros((M,1))
    
    for i in xrange(0,M):
        Z = random.gauss(0,1)
        sn = s0*numpy.exp((r-0.5*pow(v,2))*T + v*numpy.sqrt(T)*Z)
        
        sols[i] = max(0,sn-k)*numpy.exp(-r*T)
        
    
    return numpy.mean(sols), numpy.sqrt(numpy.var(sols))/numpy.sqrt(M)
    
    
def computeDelta(payoff = 'EuropeanC',k=99.0, s0=100.0, r=0.06, v=0.2, T=1, M=50000, bump = 0.1, sameSeed = True):
    
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
    
def computeUncertanty(function, M=20):
    out = numpy.zeros((M,1))
    out = [function() for x in out]
    return numpy.mean(out),numpy.var(out)
    
    
def computeDeltaApprox(k=99.0, s0=100.0, r=0.06, v=0.2, T=1, M=50000):
    s = numpy.zeros((M,1))    
    
    for i in xrange(0,M):
        Z = random.gauss(0,1)
        s[i] = s0*numpy.exp((r-0.5*pow(v,2))*T + v*numpy.sqrt(T)*Z)
    
    
    sigmoid = lambda x,k: 1/(1+numpy.exp(-(x-k)))
    sigmoid_d = lambda x,k: sigmoid(x,k)*(1-sigmoid(x,k))
    
    sigmoid_d2 = lambda x,k: numpy.abs(x-k)<0.5
    
    deltas = [sigmoid_d(x,k)*x/s0 for x in s]    
    
    return numpy.mean(deltas)
        

def computeDeltaLikelihood(k=99.0, s0=100.0, r=0.06, v=0.2, T=1, M=50000):
    shocks = numpy.zeros((M,1))
    shocks = [random.gauss(0,1) for x in shocks]
    
    
    s_t = lambda Z: s0*numpy.exp((r-0.5*pow(v,2))*T + v*numpy.sqrt(T)*Z)
    #weird_function = lambda Z: numpy.exp(-r*T)*(s_t(Z)>k)*Z/(v*s0*numpy.sqrt(T))
    weird_function = lambda Z: (max(s_t(Z)-k,0))*Z/(v*s0*numpy.sqrt(T))
    #weird_function = lambda Z: (s_t(Z)>k)*Z/(v*s0*numpy.sqrt(T))
    
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
    return mean(payoffs), numpy.sqrt(numpy.var(payoffs))/numpy.sqrt(M)
    
    
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

def plotCovarianceOfAntithetic(k=99.0, s0=100.0, r=0.06, v=0.2, T=1, M=100):
    k = numpy.linspace(10,10000,50)
    s_t = lambda Z: s0*numpy.exp((r-0.5*pow(v,2))*T + v*numpy.sqrt(T)*Z)
    
    for i in k:
        valsplus = numpy.zeros((i,1))
        #print valsplus
        valsplus = [random.gauss(0,1) for x in valsplus]
        valsminus = (-1)*numpy.array(valsplus)
        #print valsminus, valsplus
        valsplus = [s_t(z) for z in valsplus]
        valsminus = [s_t(z) for z in valsminus]
        print numpy.cov(valsplus,valsminus)
    
                

def computeAsian(k=99.0, s0=100.0, r=0.06, v=0.2, T=1, M=10000, N = 10):
    delta_t = 1.0/365
    time_steps = T/delta_t
    
    values = numpy.zeros((M,time_steps))
    values[:,0] = s0;
    for i in xrange(1,int(time_steps)):
        print i
        shocks = scipy.randn(1,M)
        shocks = shocks*v*numpy.sqrt(delta_t)+r*delta_t
        values[:,i] = (numpy.multiply(shocks,values[:,i-1])+values[:,i-1])
        print numpy.var(shocks)
        
    relevant_values = values[:,time_steps-N:]
    relevant_values = numpy.mean(relevant_values,axis=1)
    payoff = [max(x-k,0)*numpy.exp(-r*T) for x in relevant_values]
    #numpy.mean(relevant_values-k)
        #= numpy.multiply(shocks,values[:,i-1])+values[:,i-1]
    return numpy.mean(payoff)
    
    
def computeAsianCV(k=99.0, s0=100.0, r=0.06, v=0.2, T=1, M=10000, N = 10):
    delta_t = 1.0/365
    time_steps = T/delta_t
    
    values = numpy.zeros((M,time_steps))
    values[:,0] = s0;
    for i in xrange(1,int(time_steps)):
        print i
        shocks = scipy.randn(1,M)
        shocks2 = shocks*v*numpy.sqrt(delta_t)+r*delta_t
        values[:,i] = (numpy.multiply(shocks2,values[:,i-1])+values[:,i-1])
        print values[:,i].shape
        
    relevant_values = values[:,time_steps-N:]
    relevant_values_CA = numpy.mean(relevant_values,axis=1)
    relevant_values_CB = scipy.stats.gmean(relevant_values, axis=1)
    print relevant_values_CB
    print relevant_values_CA
    
    #test_N = N*(T+0.0)/365
    weird_sig = numpy.sqrt(T*pow(v,2)*(N+1)*(2*N+1)/(6*pow(N,2)))
    weird_mean = T*(r - pow(v,2)/2)*(N+1)/(2*N)
    #perfect_values_CB = scipy.randn(1,M)
    perfect_values_CB = shocks*weird_sig + weird_mean
    perfect_values_CB = numpy.exp(perfect_values_CB)
    perfect_values_CB = perfect_values_CB*s0
    print perfect_values_CB
    
    d1 = (numpy.log(s0/k)+weird_mean)/(weird_sig)
    d2 = d1 - weird_sig
    stock = s0*scipy.stats.norm.cdf(d1)
    opt = k * numpy.exp(-r*T)*scipy.stats.norm.cdf(d2)
    payoff_CA = [max(x-k,0)*numpy.exp(-r*T) for x in relevant_values_CA]
    #numpy.mean(relevant_values-k)
        #= numpy.multiply(shocks,values[:,i-1])+values[:,i-1]
    return stock-opt, numpy.mean(payoff_CA)