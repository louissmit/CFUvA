# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 11:25:55 2014

@author: henning
"""
import matplotlib.pyplot as plt
import numpy
import pprint
import inspect
from blackscholes import BS

def FD(I=100,N=500, r=0.04, v=0.3, s0 = 100.0, k = 110.0, T=1.0, M1=-2.0, M2=2.0, type='ftcs'):
	V = numpy.zeros((I+2,N+1))

	l = numpy.log(s0)

	print 'smin', numpy.exp(l+M1), ' smax',numpy.exp(l+M2)

	V[:-1,0] = numpy.linspace(M1+l, l+M2, I+1)
	delta_x = V[2,0] - V[1,0]

	print delta_x

	for i in xrange(0,I+1):
		V[i,0]=max(numpy.exp(V[i,0])-k,0)


	delta_t = T/N


	alpha = (r-pow(v,2)/2)*(delta_t/(2*delta_x))
	beta = (pow(v,2)/2)*(delta_t/pow(delta_x,2))
	gamma = r*delta_t

	if type == 'ftcs':
		a1=0
		a0=1
		a_1=0
	else:
		a1 = -alpha - beta/2
		a0 = beta + gamma/2 + 1
		a_1 = alpha - beta/2

	A = numpy.zeros((I, I))
	#FTCS:
	for i in xrange(0,I):
		A[i,i]=a0
		if i < I-1:
			A[i+1,i]=a1
			A[i,i+1]=a_1

	Abar = numpy.linalg.inv(A)




	if type == 'ftcs':
		b1 = (alpha + beta)
		b0 = (1-2*beta-gamma)
		b_1 = (beta - alpha)
	else:
		b1 = -a1
		b0 = -beta - gamma/2 + 1
		b_1 = -a_1

	print 'a', a1, a0, a_1
	print 'b', b1, b0, b_1

	for n in xrange(0,N):
		#calculate c
		c = numpy.zeros((I,1))
		#not sure if c[0] and c[I] make sense
		c[0] = b1*V[2,n] + b0*V[1,n] + b_1*V[0, n] - a_1*V[0,n+1]
		c[I-1] = b1*V[I+1,n] + b0*V[I,n]+ b_1*V[I-1,n] - a1*V[I+1, n+1]
		for j in xrange(1,I-1):
			c[j] = b1*V[j+1,n] + b0*V[j,n] + b_1*V[j-1,n]

		V[0:I ,n+1] = numpy.dot(Abar,c).T

	return V

V = FD(type='ftcs')

def plotDelta():
	args = inspect.getargspec(FD)
	N = args.defaults[args.args.index('N')]
	s0 = args.defaults[args.args.index('s0')]
	I = args.defaults[args.args.index('I')]
	M1 = args.defaults[args.args.index('M1')]
	M2 = args.defaults[args.args.index('M2')]
	VN = FD(type='ftcs')[:,N]

	l = numpy.log(s0)
	# M = numpy.gradient(VN, (M2 - M1)/(I+1))
	log_range = numpy.linspace(l+M1, l+M2, I+1)

	M = []
	for i in xrange(1, len(VN)-1):
		delta = (VN[i] - VN[i-1]) / (numpy.exp(log_range[i]) - numpy.exp(log_range[i-1]))
		M.append(delta)

	S = []
	for i in xrange(1,len(log_range)-1):
		S.append((numpy.exp(log_range[i]) + numpy.exp(log_range[i-1])) / 2)
	B = [BS(s0=s)[1] for s in S]
	plt.plot(M, label="FTCS")
	plt.plot(B, label="Black-Scholes")
	plt.legend()
	plt.show()

	return M, B, VN

M, B, VN = plotDelta()

