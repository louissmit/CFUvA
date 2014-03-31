# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 14:50:42 2014

@author: henning
"""

import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import cmath
from blackscholes import BS
from ex3 import FD
from matplotlib import cm
from matplotlib import colors

#not sure if t = T
phi = lambda u, r, sig, T: np.exp(1j * u * (r - pow(sig, 2) / 2) * T - pow(sig, 2) * T * pow(u, 2) / 2)


def F(n, a, b, r, sig, T):
	u = n * np.pi / (b - a)
	res = 2 / (b - a) * np.real(phi(u, r, sig, T) * np.exp(-1j * a * u))
	if n == 0:
		return res / 2
	return res


def G(n, a, b, K):
	weirdX = (1 / (1 + pow((n * np.pi) / (b - a), 2))) * (np.cos(n * np.pi) * np.exp(b) -
														  np.cos((-n * np.pi * a / (b - a))) +
														  n * np.pi / (b - a) * np.sin(n * np.pi) * np.exp(b) -
														  n * np.pi / (b - a) * np.sin(-n * np.pi * a / (b - a)))

	# print "chi " + str(weirdX)
	#weirdV = 0
	if n == 0:
		weirdV = b
	else:
		weirdV = (b - a) / (n * np.pi) * (np.sin(n * np.pi) - np.sin(-n * np.pi * a / (b - a)))

	return K * (weirdX - weirdV) * 2 / (b - a)


def cos_shit(n=124, s0=100.0, K=100.0, T=1.0, sig=0.4, r=0.03):
	a = np.log(s0 / K) + r * T - 12 * np.sqrt(pow(sig, 2) * T)
	b = np.log(s0 / K) + r * T + 12 * np.sqrt(pow(sig, 2) * T)

	X = np.log(s0 / K)
	summ = 0
	for i in xrange(0, n):
		# print str(G(i)) +"\t"+ str(F(i))
		summ += G(i, a, b, K) * F(i, a, b, r, sig, T)
	return np.exp(-r * T) * (b - a) / 2 * summ


def fourier_coefficients():
	K = xrange(1, 7)
	b = 5.0
	a = -5.0
	phi = lambda u: np.exp(-pow(u, 2) / 2.0)
	error = []
	X = np.arange(-5, 5, 0.1)
	print X
	gold = norm.pdf(X)
	F = lambda n, a, b: ( 2 / (b - a)) * np.real(
		phi(n * np.pi / (b - a)) * np.exp(-cmath.sqrt(-1) * (n * a * np.pi / (b - a))))
	for k in K:
		N = pow(2, k)
		f = 0.5 * F(0, a, b)
		for n in xrange(1, N):
			f += F(n, a, b) * np.cos(n * np.pi * (X - a) / (b - a))

		plt.plot(X, f, label="N = " + str(N))
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
	plt.xlim([0, 7])
	plt.show()


def plot_convergence(n=124, s0=100.0, K=100.0, T=1.0, sig=0.4, r=0.03):
	Ns = [16, 32, 64, 96, 128, 160, 1920]
	Ns = xrange(6, 120, 1)
	Y = []
	error = []
	for N in Ns:
		value = cos_shit(N)
		# Y.append(value)
		error.append(value - BS(k=K, s0=s0, r=r, v=sig)[0])

	# BY = [BS(k=K,s0=s0,r=r,v=sig)[0] for i in xrange(0, len(Y))]

	plt.plot(Ns, error, label='Error w.r.t Black-Scholes')
	plt.xlabel('$N$ Fourier-Cosine coefficients')
	plt.ylabel('Error')
	plt.xlim([6, 120])
	plt.legend()
	plt.show()


def calculate_grid(I=100, N=100, s0=100.0, K=100.0, T=1.0, sig=0.4, r=0.03):
	Vfd, s0s = FD(I=I, N=N, v=sig, s0=s0, k=K)
	Vbs = np.zeros((I + 1, N + 1))
	Vcos = np.zeros((I + 1, N + 1))
	delta_t = T/N
	for i in xrange(0, I):
		for n in xrange(0, N):
			if n == 0:
				Vbs[i,n] = max(s0s[i]-K,0)
				Vcos[i,n] = max(s0s[i]-K,0)
			else:
				Vcos[i, n] = cos_shit(n=25, T=delta_t*n, s0=s0s[i])
				Vbs[i, n] = BS(k=K,s0=s0s[i],r=r,v=sig, T=delta_t*n)[0]
	return s0s, Vcos,Vbs,Vfd

def plotV(N, V, s0s):
	n = np.linspace(0, N+1, N+1)
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	X, Y = np.meshgrid(n, s0s[0:60])

	norm = colors.Normalize(vmin = 0, vmax = 50, clip = True)
	ax.plot_surface(X,Y, V[0:60,:], cmap=cm.coolwarm, norm=norm)
	ax.set_xlabel("$n$")
	ax.set_ylabel("$S_o$")
	ax.set_zlabel("$V$")

# fourier_coefficients()
# print cos_shit()
# plot_convergence()

s0s, V, Vbs, Vfd = calculate_grid()
plotV(100, V, s0s)
plotV(100, Vbs, s0s)
plotV(100, Vfd, s0s)
plt.show()
