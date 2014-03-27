
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import cmath

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

fourier_coefficients()
