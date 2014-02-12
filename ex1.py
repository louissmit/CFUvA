# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle
from blackscholes import BS


def computeMat(k=99, s0=100, r = 0.06, v=0.2, N=50, type="call", european=True):
	dt = 1.0 / N
	u = np.exp(v*np.sqrt(dt))
	d = np.exp(-v*np.sqrt(dt))

	p = (np.exp(r*dt)-d)/(u-d)

	vals = np.zeros((N,N))
	S = np.zeros((N,N))

	#initliaze vals
	for j in xrange(0,N):
		s = s0*pow(u,j)*pow(d,N-1-j)
		#check option type- put or call
		if type == "call":
			vals[j,N-1] = max(0, s - k)
		elif type == "put":
			vals[j,N-1] = max(0, k - s)
		else:
			vals[j,N-1] = k-s
	for i in range(N-2,-1,-1):
		for j in xrange(0,N-1):
			s = s0*pow(u,j)*pow(d,i-j)
			S[j, i] = s
			f = np.exp(-r*dt)*(p*vals[j+1,i+1] + (1-p)*vals[j,i+1])
			#check option type - european or american
			if european:
				vals[j,i] = f
			else:
				if type == "call":
					vals[j,i] = max(s - k, f)
				elif type == 'put':
					vals[j,i] = max(k - s, f)
				else:
					vals[j,i] = f

	return vals, S


def binomialConvergence():
	X = []
	Y = []
	for i in xrange(5, 150, 1):
		X.append(i)
		price = computeMat(N=i)[0,0]
		Y.append(price)
	pickle.dump(Y, open( "prices", "wb" ))
	plt.plot(X, Y)
	plt.xlabel('iterations N')
	plt.ylabel('call option price in euro')
	plt.plot(X, [BS() for x in X])
	plt.show()

def plotPriceVSVolatility():
	X = []
	Y = []
	Y2 = []
	vold1 = np.linspace(0.01, 0.15, 10)
	for v in vold1:
		X.append(v * 100)
		vals, S = computeMat(v=v, type='call', european=True)
		Y.append(vals[0,0])
		s, deltabs = BS(vd1=v, vd2=v)
		Y2.append(s)

	plt.plot(X, Y, label="Binomial tree")
	plt.plot(X, Y2, label="Black-Scholes")
	plt.xlabel('volatility in %')
	plt.ylabel('call option price in euro')
	plt.legend(loc=2)

	plt.show()

def plotDeltaVSVolatility():
	X = []
	Y = []
	Y2 = []
	vold1 = np.linspace(0.01, 0.9, 50)
	for v in vold1:
		X.append(v * 100)
		vals, S = computeMat(v=v, type='call', european=True)

		delta = (vals[1, 1] - vals[1,2]) / (S[1,1] - S[1, 2])
		Y.append(delta)
		s, deltabs = BS(vd1=v, vd2=v)
		Y2.append(deltabs)

	plt.plot(X, Y, label="Binomial tree")
	plt.plot(X, Y2, label="Black-Scholes")
	plt.xlabel('volatility in %')
	plt.ylabel('delta')
	plt.legend()

	plt.show()
# binomialConvergence()
# plotPriceVSVolatility()
plotDeltaVSVolatility()

# print computeMat(european=False, type="call")[0,0]
# print computeMat(european=False, type="put")[0,0]
# print computeMat(european=False, type="forward")[0,0]
