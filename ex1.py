# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle
from blackscholes import BS


def computeMat(k=99.0, s0=100.0, r = 0.06, v=0.2, N=50, type="call", european=True):
	dt = 1.0 / N
	u = np.exp(v*np.sqrt(dt))
	d = np.exp(-v*np.sqrt(dt))

	p = (np.exp(r*dt)-d)/(u-d)

	vals = np.zeros((N,N+1))
	S = np.zeros((N,N+1))

	#initliaze vals
	for j in xrange(0,N):
		s = s0*pow(u,j)*pow(d,N-1-j)
		#check option type- put or call
		if type == "call":
			vals[N-1, j] = max(0, s - k)
		elif type == "put":
			vals[N-1, j] = max(0, k - s)
		else:
			vals[N, j] = k-s
	for i in range(N-2, -1,-1):
		for j in xrange(0,N):
			s = s0*pow(u,j)*pow(d,i-j)
			S[i, j] = s
			f = np.exp(-r*dt)*(p*vals[i+1,j+1] + (1-p)*vals[i+1,j])
			#check option type - european or american
			if european:
				vals[i,j] = f
			else:
				if type == "call":
					vals[i,j] = max(s - k, f)
				elif type == 'put':
					vals[i,j] = max(k - s, f)
				else:
					vals[i,j] = f

	return vals, S


def binomialConvergence():
	X = []
	Y = []
	Y2 = []
	for i in xrange(1, 101, 1):
		vals, S = computeMat(N=i)
		price = vals[0,0]

		X.append(i)
		Y.append(price)
		# if i % 2 == 0:
		# 	X.append(i)
		# 	Y.append(price)
		# else:
		# 	Y2.append(price)
	# pickle.dump(Y, open( "prices", "wb" ))

	# plt.plot(X, Y, label="even")
	plt.plot(X, Y, label="binomial tree")
	# plt.plot(X, Y2, label="odd")
	plt.xlabel('N')
	plt.ylabel('call option price in euro')
	plt.plot(X, [BS()[0] for x in X], label="Black-Scholes")
	plt.legend(loc=4)
	plt.show()

def plotPriceVSVolatility():
	X = []
	Y = []
	Y2 = []
	Y3 = []
	vold1 = np.linspace(0.01, 0.9, 10)
	for v in vold1:
		X.append(v * 100)
		vals, S = computeMat(v=v, type='put', european=True)
		Y.append(vals[0,0])
		# s, deltabs = BS(vd1=v, vd2=v)
		# Y2.append(s)
		vals, S = computeMat(v=v, type='put', european=False)
		Y3.append(vals[0,0])

	plt.plot(X, Y, label="Binomial tree")
	# plt.plot(X, Y2, label="Black-Scholes")
	plt.plot(X, Y3, label="american")
	plt.xlabel('volatility in %')
	plt.ylabel('call option price in euro')
	plt.legend(loc=2)

	plt.show()

def plotDeltaVSVolatility(K=[99]):

	vold1 = np.linspace(0.01, 0.9, 50)
	i = 1
	for k in K:
		X = []
		Y = []
		Y2 = []
		for v in vold1:
			X.append(v * 100)
			vals, S = computeMat(v=v, type='call', european=True, k=k)

			delta_f = vals[1, 1] - vals[2,1]
			delta_S = S[1, 1] - S[2,1]
			delta = delta_f / delta_S
			Y.append(delta)
			s, deltabs = BS(vd1=v, vd2=v, k=k+0.0)
			Y2.append(deltabs)

		plot = plt.subplot(1 , len(K),i)
		plot.plot(X, Y, label="Binomial tree")
		plot.plot(X, Y2, label="Black-Scholes")
		plot.set_xlabel('volatility in %')
		plot.set_ylabel('delta for k=' + str(k))
		plot.legend()
		i+=1

	plt.show()

# binomialConvergence()
# plotPriceVSVolatility()
plotDeltaVSVolatility(K=[50, 99, 120, 150])

# vals,s = computeMat(european=True, type="call")
# print vals[0,0]
# vals,s = computeMat(european=False, type="call")
# print vals[0,0]
# print computeMat(european=False, type="put")[0,0]
# print computeMat(european=False, type="forward")[0,0]
