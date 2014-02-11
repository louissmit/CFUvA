# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle


def computeMat(k=99, s0=100, r = 0.06, v=0.2, N=100):
	dt = 1.0 / N
	u = np.exp(v*np.sqrt(dt))
	d = np.exp(-v*np.sqrt(dt))

	p = (np.exp(r*dt)-d)/(u-d)

	vals = np.zeros((N,N))

	#initliaze vals
	for j in xrange(0,N):
		s = s0*pow(u,j)*pow(d,N-1-j)
		vals[j,N-1] = max(0,s-k)
	for i in range(N-2,-1,-1):
		for j in xrange(0,N-1):
			vals[j,i] = np.exp(-r*dt)*(p*vals[j+1,i+1] + (1-p)*vals[j,i+1])

	return vals

X = []
for i in xrange(5, 1000, 10):
	price = computeMat(N=i)[0,0]
	print price
	X.append(price)


pickle.dump(X, open( "prices", "wb" ))
plt.plot(X)
plt.show()
