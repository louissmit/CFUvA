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
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import colors

def FD(I=100,N=500, r=0.06, v=0.2, s0 = 100.0, k = 99.0, T=1.0, M1=-2.0, M2=2.0, type='ftcs'):
	V = numpy.zeros((I+1,N+1))

	l = numpy.log(s0)
	V[:,0] = numpy.linspace(M1+l, l+M2, I+1)

	s0s = numpy.exp(numpy.linspace(M1+l, l+M2, I+1))

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
	elif type == 'cn':
		a1 = -alpha - beta/2
		a0 = beta + gamma/2 + 1
		a_1 = alpha - beta/2
	else:
		raise TypeError("not a valid scheme type")

	A = numpy.zeros((I+1, I+1))
	#FTCS:
	for i in xrange(0,I+1):
		A[i,i]=a0
		if i < I:
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
		c = numpy.zeros((I+1,1))
		#not sure if c[0] and c[I] make sense
		if type=='ftcs':
			c[0] = b1*V[2,n] + b0*V[1,n] + b_1*V[0, n] - a_1*V[0,n+1]
			c[I] = b0*V[I,n]+ b_1*V[I-1,n] + b1*(V[I,n]+delta_x*numpy.exp(l+M2))#- a1*V[I+1, n+1] + b1*V[I+1,n]
		else:
			c[0] = b1*V[2,n] + b0*V[1,n] + b_1*V[0, n] - a_1*V[0,n+1]
			# c[I] = b0*V[I,n]+ beta*V[I-1,n] + b1*(V[I,n]+delta_x*numpy.exp(l+M2))
			c[I] = (b1-a1)*(numpy.exp(-r*(T-n*delta_t))*max(numpy.exp(l+M2)-k,0)) + b0*V[I,n] + b_1*V[I-1,n]

		for j in xrange(1,I):
			c[j] = b1*V[j+1,n] + b0*V[j,n] + b_1*V[j-1,n]

		V[:,n+1] = numpy.dot(Abar,c).T

	return V, s0s


#V, s = FD(type='ftcs')

def plotDelta(**kwargs):
	args = inspect.getargspec(FD)
	N = kwargs.get('N') or args.defaults[args.args.index('N')]
	s0 = kwargs.get('s0') or args.defaults[args.args.index('s0')]
	I = args.defaults[args.args.index('I')]
	M1 = args.defaults[args.args.index('M1')]
	M2 = args.defaults[args.args.index('M2')]
	kwargs['type'] = 'ftcs'
	V_ftcs = FD(**kwargs)[0][:,N]
	kwargs['type'] = 'cn'
	V_crank = FD(**kwargs)[0][:,N]

	l = numpy.log(s0)
	# M = numpy.gradient(VN, (M2 - M1)/(I+1))
	log_range = numpy.linspace(l+M1, l+M2, I+1)

	M = []
	C = []
	for i in xrange(1, len(V_ftcs)-1):
		delta = (V_ftcs[i] - V_ftcs[i-1]) / (numpy.exp(log_range[i]) - numpy.exp(log_range[i-1]))
		delta_crank = (V_crank[i] - V_crank[i-1]) / (numpy.exp(log_range[i]) - numpy.exp(log_range[i-1]))
		M.append(delta)
		C.append(delta_crank)

	B = []
	S = []
	kwargs.pop('type')
	for i in xrange(1,len(log_range)-1):
		s = (numpy.exp(log_range[i]) + numpy.exp(log_range[i-1])) / 2
		kwargs['s0'] = s
		S.append(s)
		B.append(BS(**kwargs)[1])

	plt.plot(S, C, label="Crank-Nicholson")
	plt.plot(S, M, label="FTCS")
	plt.plot(S, B, label="Black-Scholes")
	plt.xlabel("$S_0$")
	plt.ylabel("$\delta$")
	plt.legend(loc=4)
	plt.show()

	return M, B, V_ftcs

# M, B, VN = plotDelta(k=110.0, s0=110.0, v=0.3, r=0.04)


def plotV(N1=1000,I1=100):
	mat,s = FD(N=N1, I=I1, type='cn', k=110.0, s0=110.0, v=0.3, r=0.04)
	cn_mat,s = FD(N=N1, I=I1, type='cn')
	x = [BS(s0=s0)[0] for s0 in s]

	# s = s[0:-1]
	# x = x[0:-1]

	# y = mat[0:-1,-1]
	# y_cn = cn_mat[0:-1,-1]


	n = numpy.linspace(0, N1+1, N1+1)
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	X, Y = numpy.meshgrid(n,s[0:60])
	print numpy.shape(X),numpy.shape(Y),numpy.shape(mat)

	norm = colors.Normalize(vmin = 0, vmax = 50, clip = True)
	ax.plot_surface(X,Y, mat[0:60,:], cmap=cm.coolwarm, norm=norm)
	ax.set_xlabel("$n$")
	ax.set_ylabel("$S_o$")
	ax.set_zlabel("$V$")

	# plt.plot(s,x,color='blue', label="Black-Scholes")
	# plt.plot(s,y,color='green', label="FTCS")
	# plt.plot(s,y_cn,color='red', label="Crank-Nicholson")
	# plt.xlabel("$S_0$")
	# plt.ylabel("$V$")
	# plt.legend(loc=4)
	plt.show()

# M,s = FD()
# plotV()