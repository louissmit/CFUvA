# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle


def computeMat(k=99, s0=100, r = 0.06, v=0.2, N=50, type="call", european=True):
	dt = 1.0 / N
	u = np.exp(v*np.sqrt(dt))
	d = np.exp(-v*np.sqrt(dt))

	p = (np.exp(r*dt)-d)/(u-d)

	vals = np.zeros((N,N))

	#initliaze vals
	for j in xrange(0,N):
		s = s0*pow(u,j)*pow(d,N-1-j)
		#check option type- put or call
		if type == "call":
			vals[j,N-1] = max(0, s - k)
		elif type == "put":
			vals[j,N-1] = max(0, k-s)
		else:
			vals[j,N-1] = k-s
	for i in range(N-2,-1,-1):
		for j in xrange(0,N-1):
			#check option type - european or american
			if european:
				vals[j,i] = np.exp(-r*dt)*(p*vals[j+1,i+1] + (1-p)*vals[j,i+1])
			else:
				s = s0*pow(u,j)*pow(d,i-j)
				f = np.exp(-r*dt)*(p*vals[j+1,i+1] + (1-p)*vals[j,i+1])
				print s - k
				print f
				print "----"
				if type == "call":
					vals[j,i] = max(s - k, f)
				elif type == 'put':
					vals[j,i] = max(k - s, f)
				else:
					vals[j,i] = f

	return vals


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
	plt.show()

def plotPriceVSVolatility(call=True, european=True):
	X = []
	Y = []
	Y2 = []
	vold1 = np.linspace(0.01, 0.15, 10)
	for v in vold1:
		X.append(v * 100)
		Y.append(computeMat(v=v, call=call, european=european)[0,0])
		Y2.append(computeMat(v=v, call=call)[0,0])

	plt.plot(X, Y)
	plt.plot(X, Y2)
	plt.xlabel('volatility in %')
	plt.ylabel('call option price in euro')

	plt.show()

# binomialConvergence()
# plotPriceVSVolatility(european=False)

print computeMat(european=False, type="call")[0,0]
print computeMat(european=False, type="put")[0,0]
print computeMat(european=False, type="forward")[0,0]
