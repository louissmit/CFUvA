# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 14:50:42 2014

@author: henning
"""

import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import cmath

s0 = 100.0
K = 99.0
T = 1.0
sig = 0.2
r = 0.06
a = np.log(s0/K) + r*T - 12*np.sqrt(pow(sig,2)*T)
b = np.log(s0/K) + r*T + 12*np.sqrt(pow(sig,2)*T)

X = np.log(s0/K)

#not sure if t = T
phi = lambda u: np.exp(1j*u*(r-pow(sig,2)/2)*T - pow(sig,2)*T*pow(u,2)/2)

def F(n):
# F = lambda n, a, b : ( 2 / (b-a)) * np.real(phi(n*np.pi / (b - a))*np.exp(-cmath.sqrt(-1)*(n*a*np.pi / (b - a))))
	return 2/(b-a)*np.real(phi(n*np.pi/(b-a))*np.exp(-1j*n*a*np.pi/(b-a)))

def G(n):
	weirdX = 1/(1+pow((n*np.pi)/(b-a),2))*(np.cos(n*np.pi)*np.exp(b) -
										   np.cos((n*np.pi*a/(b-a))) + n*np.pi/(b-a)*np.sin(n*np.pi)*np.exp(b)
										   - n*np.pi/(b-a)*np.sin(n*np.pi*a/(b-a)))

	#weirdV = 0
	if n == 0:
		weirdV = b
	else:
		weirdV = (b-a)/(n*np.pi)*(np.sin(n*np.pi) - np.sin(n*np.pi*a/(b-a)))

	return K*(weirdX-weirdV)*2/(b-a)

def cos_shit(n=124):
	summ = 0
	for i in xrange(0,n):
		summ += G(i)*F(i)
	return np.exp(-r*T)*(b-a)/2*summ




def fourier_coefficients():
	K = xrange(1,7)
	b = 5.0
	a = -5.0
	phi = lambda u: np.exp(-pow(u,2)/2.0)
	error = []
	X = np.arange(-5, 5 ,0.1)
	print X
	gold = norm.pdf(X)
	F = lambda n, a, b : ( 2 / (b-a)) * np.real(phi(n*np.pi / (b - a))*np.exp(-cmath.sqrt(-1)*(n*a*np.pi / (b - a))))
	for k in K:
		N = pow(2,k)
		f = 0.5 * F(0,a,b)
		for n in xrange(1,N):
			f += F(n, a, b) * np.cos(n*np.pi*(X - a)/(b - a))

		plt.plot(X, f, label="N = "+str(N))
		error.append(np.linalg.norm(abs(gold - f)))
	plt.plot(X, gold, label="standard normal")
	plt.legend()
	plt.xlim([-5, 5])
	plt.show()
	plt.plot(K, error)
	plt.xlabel("N")
	plt.ylabel("error")
	plt.xscale('log')
	plt.yscale('log')
	plt.xlim([0,7])
	plt.show()

# fourier_coefficients()
print cos_shit()