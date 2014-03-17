from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt

def BS(k=99.0, s0=100.0, r = 0.06, v=0.2, T=1, option_type='call'):
	# d1 = (np.log(s0/k)+(r+pow(vd1,2))*T)/(vd1*np.sqrt(T))
	d1 = (np.log(s0/k)+(r+0.5*pow(v,2))*T)/(v*np.sqrt(T))
	d2 = d1 - v/np.sqrt(T)
	stock = s0*norm.cdf(d1)
	opt = k * np.exp(-r*T)*norm.cdf(d2)
	opt_price = 0
	if option_type == 'call':
		opt_price = stock - opt
	elif option_type == 'digital':
		opt_price = np.exp(-r*T)*norm.cdf(d2)
	else:
		raise TypeError("not a valid option type")

	return opt_price, norm.cdf(d1), norm.cdf(d2)

# print BS(type='digital')
def plotPriceVSVolatility():
	X = []
	Y = []
	Y2 = []
	Y3 = []
	vold1 = np.linspace(0.01, 0.9, 30)
	for v in vold1:
		X.append(v * 100)
		s, cdf1, cdf2 = BS(vd1=v, vd2=v)
		Y.append(s)
		Y2.append(cdf1)
		Y3.append(cdf2)
		# vals, S = computeMat(v=v, type='put', european=False)
		# Y3.append(vals[0,0])

	# plt.plot(X, Y, label="Black-Scholes")
	plt.plot(X, Y2, label="cdf d1")
	plt.plot(X, Y3, label="cdf d2")
	plt.xlabel('volatility in %')
	plt.ylabel('call option price in euro')
	plt.legend(loc=2)

	plt.show()

# plotPriceVSVolatility()