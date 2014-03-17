# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/henning/.spyder2/.temp.py
"""

import numpy as np
from blackscholes import BS
import random
import matplotlib.pyplot as plt


def computeMCValuation(k=99.0, s0=100.0, r=0.06, v=0.2, T=1, M=100000, type="call"):
	sols = np.zeros((M,1))

	for i in xrange(0,M):
		Z = random.gauss(0,1)
		sn = s0*np.exp((r-0.5*pow(v,2))*T + v*np.sqrt(T)*Z)

		if type=="call":
			sols[i] = max(0,sn-k)*np.exp(-r*T)
		elif type=="put":
			sols[i] = max(0,k-sn)*np.exp(-r*T)
		else:
			raise TypeError("not a valid option type")

	return np.mean(sols), np.sqrt(np.var(sols))/np.sqrt(M)

def getPrice(**kwargs):
	return computeMCValuation(**kwargs)[0]

def computeDelta(type = 'EuropeanC',k=99.0, s0=100.0, r=0.06, v=0.2, T=1, M=100000, bump = 0.01, sameSeed = True):

	shocks = np.zeros((M,1))

	for i in xrange(0,M):
		Z = random.gauss(0,1)
		shocks[i] = np.exp((r-0.5*pow(v,2))*T + v*np.sqrt(T)*Z)

	s = np.array(shocks)*s0

	if not sameSeed:
		for i in xrange(0,M):
			Z = random.gauss(0,1)
			shocks[i] = np.exp((r-0.5*pow(v,2))*T + v*np.sqrt(T)*Z)


	s_bumped = np.array(shocks)*(s0+bump)

	#computing payoffs
	if type == 'EuropeanC':
		s_payoff = [max(x-k,0) for x in s]
		sb_payoff = [max(x-k,0) for x in s_bumped]
	elif type == 'Digital':
		s_payoff = [x>k for x in s]
		sb_payoff = [x>k for x in s_bumped]


	return (np.mean(sb_payoff)-np.mean(s_payoff))/bump


def computeDeltaApprox(k=99.0, s0=100.0, r=0.06, v=0.2, T=1, M=100000):
	s = np.zeros((M,1))

	for i in xrange(0,M):
		Z = random.gauss(0,1)
		s[i] = s0*np.exp((r-0.5*pow(v,2))*T + v*np.sqrt(T)*Z)


	sigmoid = lambda x,k: 1/(1+np.exp(-(x-k)))
	sigmoid_d = lambda x,k: sigmoid(x,k)*(1-sigmoid(x,k))

	# sigmoid_d2 = lambda x,k: np.abs(x-k)<0.5

	deltas = [sigmoid_d(x,k)*x/s0 for x in s]

	return np.mean(deltas)


def computeDeltaLikelihood(k=99.0, s0=100.0, r=0.06, v=0.2, T=1, M=100000):
	shocks = np.zeros((M,1))
	shocks = [random.gauss(0,1) for x in shocks]


	s_t = lambda Z: s0*np.exp((r-0.5*pow(v,2))*T + v*np.sqrt(T)*Z)
	weird_function = lambda Z: np.exp(-r*T)*(s_t(Z)>k)*Z/(v*s0*np.sqrt(T))

	deltas = [weird_function(Z) for Z in shocks]

	return np.mean(deltas)


def plotValuationConvergence(n=10, ms = np.linspace(50,100000,15)):
	X=[]
	Y=[]

	for m in ms:
		for j in xrange(0,n):
			#print int(m)
			X.append(int(m))
			Y.append(computeMCValuation(M=int(m))[0])

	plt.scatter(X,Y)

def varianceAntithetic(k=99.0, s0=100.0, r=0.06, v=0.2, T=1, M=100, meanPayoff=True):
	shocks = np.zeros((M,1))
	shocks = [random.gauss(0,1) for x in shocks]
	#shocks = [max(x,-x) for x in shocks]
	shocks2= np.array(shocks)*(-1)

	#shocks = np.concatenate((shocks,shocks2),axis=1)


	s_t = lambda Z: s0*np.exp((r-0.5*pow(v,2))*T + v*np.sqrt(T)*Z)

	sn = [s_t(Z) for Z in shocks]
	sn2 = [s_t(Z) for Z in shocks2]

	if not meanPayoff:
		sn = (np.array(sn)+np.array(sn2))/2

	if meanPayoff:
		payoffs = [max(0,x-k)*np.exp(-r*T) for x in sn]
		payoffs2 = [max(0,x-k)*np.exp(-r*T) for x in sn2]

		payoffs = (np.array(payoffs)+np.array(payoffs2))/2
	else:
		payoffs = [max(0,x-k)*np.exp(-r*T) for x in sn]

	#return sn
	return np.mean(payoffs), np.sqrt(np.var(payoffs))/np.sqrt(M)


def plotAntitheticShit():
	total_M = 10000
	m = np.linspace(50,5000,50)

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
		Y.append(np.mean(results1))
		Y2.append(np.mean(results2))
		Y3.append(np.mean(results3))
		X.append(mi)

	plt.plot(X,Y,label='Antithetic Payoffmean')
	plt.plot(X,Y,label='Antithetic Stockmean')
	plt.plot(X,Y2,label='Vanilla')
	plt.legend()
	plt.show()

def plotStandardError(maxM=50000):
	X = []
	Y = []
	for M in xrange(10, maxM, 100):
		X.append(M)
		std_error = computeMCValuation(M=M)[1]
		Y.append(std_error)

	plt.plot(X, Y)
	plt.legend()
	plt.show()

def plotConvergence(maxM=50000):
	X = []
	Y = []
	for M in xrange(10, maxM, 100):
		X.append(M)
		value = computeMCValuation(M=M)[0]
		Y.append(value)

	plt.plot(X, [BS()[0] for x in X], label="Black-Scholes")
	plt.plot(X, Y)
	plt.xlabel('M')
	plt.ylabel('price')
	plt.legend()
	plt.show()

def plotPriceVSVolatility():
	X = []
	Y = []
	Y2 = []
	vold1 = np.linspace(0.01, 0.9, 30)
	for v in vold1:
		X.append(v * 100)
		c = computeMCValuation(v=v)[0]
		Y.append(c)
		s = BS(vd1=v, vd2=v)[0]
		Y2.append(s)

	plt.plot(X, Y, label="Monte Carlo")
	plt.plot(X, Y2, label="Black-Scholes")
	# plt.plot(X, Y3, label="american")
	plt.xlabel('volatility in %')
	plt.ylabel('call option price in euro')
	plt.legend(loc=2)

	plt.show()


def priceVSStrike():
	P = []
	C = []
	X = []
	for k in xrange(50,150, 5):
		X.append(k)
		P.append(getPrice(type="put",k=k))
		C.append(getPrice(type="call",k=k))
	plt.plot(X, P, label="put")
	plt.plot(X, C, label="call")
	plt.xlabel('strike')
	plt.ylabel('option price')
	plt.legend()
	plt.show()

def deltaVSEpsilonM():

	for m in [pow(10,4), pow(10,5)]:
		s = ""
		for e in [0.01, 0.02, 0.5]:
			error = np.round(abs(1 - (BS()[1] / computeDelta(M=m, bump=e, sameSeed=True)))*100,2)
			s += str(error) + "\% & "

		print s + "\n"

def plotDelta(maxM = 100000):
	X = []
	Y = []
	for M in xrange(1000, maxM, 2000):
		X.append(M)
		delta = computeDelta(M=M)
		Y.append(delta)
	plt.plot(X, Y)
	plt.xlabel('M')
	plt.ylabel('delta')
	plt.show()


def deltaUncertainty(approx_type="bump", M=50000, **kwargs):
	n = 300
	X = np.empty((n,1))
	for i in xrange(0, n):
		if approx_type == 'bump':
			delta = computeDelta(M=M, **kwargs)
		elif approx_type == 'smooth':
			delta = computeDeltaApprox(M=M, **kwargs)
		elif approx_type == 'likelihood':
			delta = computeDeltaLikelihood(M=M)

		X[i] = delta
	print np.mean(X)
	print np.sqrt(np.var(X))

# deltaUncertainty(approx_type='smooth')
# print computeMCValuation(type="put")
# plotAntitheticShit()
# priceVSStrike()
# print computeDelta(M=200000)
# print BS()[1]
# print computeDelta(sameSeed=False)
# deltaVSEpsilonM()
# plotDelta()
# plotStandardError()
# plotConvergence()
plotPriceVSVolatility()
