from scipy.stats import norm
import numpy as np

def BS(k=99.0, s0=100.0, r = 0.06, vd1=0.2, vd2=0.2, T=1):
    d1 = (np.log(s0/k)+(r+pow(vd1,2))*T)/(vd1*np.sqrt(T))
    d2 = (np.log(s0/k)+(r+pow(vd2,2))*T)/(vd2*np.sqrt(T)) - vd2/np.sqrt(T)
    stock = s0*norm.cdf(d1)
    opt = k * np.exp(-r*T)*norm.cdf(d2)
    return stock-opt, norm.cdf(d1)

print BS()