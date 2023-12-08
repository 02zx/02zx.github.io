import numpy as np
data=np.loadtxt('energy.dat')
T=443
N=360
deltaA=(data[0]*0.05+data[-1]*0.05 + np.sum(data[1:-1]*0.1))*1000/N/8.314/T
print(deltaA)
