# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/henning/.spyder2/.temp.py
"""

import numpy as np
import matplotlib.pyplot as plt
import random


def computeMCValuation(k=99, s0=100.0, r=0.06, v=0.2, T=1, M=100000):
	sols = np.zeros((M,1))
	for i in xrange(0,M):
		Z = random.gauss(0,1)
		sn = s0*np.exp((r-0.5*pow(v,2))*T + v*np.sqrt(T)*Z)
		sols[i] = max(0,sn-k)*np.exp(-r*T)
	return np.mean(sols), np.sqrt(np.var(sols))/np.sqrt(M)

def plotStandardError(maxM=20000):
	X = []
	Y = []
	for M in xrange(10, maxM):
		X.append(M)
		std_error = computeMCValuation(M=M)[1]
		Y.append(std_error)

	plt.plot(X, Y)
	plt.show()

def plotConvergence(maxM=20000):
	X = []
	Y = []
	for M in xrange(10, maxM, 10):
		X.append(M)
		value = computeMCValuation(M=M)[0]
		Y.append(value)

	plt.plot(X, Y)
	plt.show()

# plotStandardError()
plotConvergence()
# print computeMCValuation()