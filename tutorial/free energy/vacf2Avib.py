#import matplotlib.pyplot as plt
import numpy as np
import sys

N_mol=288  # number of molecules sys.argv[1]
Temp=10 #simulation temperature
OW=np.loadtxt('{}/O.xvg'.format(N_mol))
HW=np.loadtxt('{}/HW.xvg'.format(N_mol))
mass_OW=15.9996
mass_HW=1.0008
period=OW[-1,0]
mwvacf=np.round(N_mol*mass_OW*OW[:,1]+2*N_mol*mass_HW*HW[:,1],5) #compute mass-weighted vacf
with open('{}/water.mwvacf'.format(N_mol),'w') as f:
    for i in range(len(OW[:,1])):
        print(OW[i,0],mwvacf[i],file=f)


def vdos(N_mol):
    filename = "./{}/water.mwvacf".format(N_mol)
    data = np.loadtxt(filename)
    
    x = data[:, 0]
    y = data[:, 1]
    
    
    N=len(y)
    print(N_mol)
    dt = 0.01
    L = dt*float(N_mol-1)
    
    print(L)
    
    for i in range(len(y)):
        if i == 0:
            continue
        y[i] = y[i] * np.sin(np.pi * i / (N_mol-1)) / (np.pi * i / (N_mol-1)) 
    
    n_values = []
    an_values = []
    
    n = 0
    while n < 500:
        i = 0
        an = 0
        while i < len(y):
    
            dx=dt
            ## A_0 = (1/T)*sum{f(t)cos(0)dt} ; cox(0)=1
            ## A_n =(2/T)sum{f(t) cos(2pi t /T)dt}
            if(i==0):
                dx=dt/2.0 #sum{f(t)cos(0)dt}/2 
            an += y[i] * np.cos(x[i]*n*2*np.pi/L)*dx #sum{f(t) cos(2pi t /T)dt}
            
            i += 1
    
        an = an *2.0 / L #(2/T)sum{f(t) cos(2pi t /T)dt}
        if n == 0:
            an /= 2.0
        
        n_values.append(n)
        an_values.append(an)
        n += 1
        
    return np.stack((np.array(n_values),np.array(an_values)),axis=0)

out_vdos=vdos(N_mol)
total=np.sum(out_vdos[1,:])*out_vdos[0,1] #integral the vdos
with open('{}/water.vdos'.format(N_mol),'w') as f:
    for i in range(out_vdos.shape[1]):
        print('#frequency(cm^-1) normalized_vdos, n, vdos',file=f)
        print(out_vdos[0,i]/period/(3/100),out_vdos[1,i]/total,out_vdos[1,i],file=f)


#degree of freedom
N_f = 6*N_mol-6

def F(freq, num,T,dv):
    return np.sum(8.314*T*num*np.log(2*np.sinh(6.626*6.02*10*freq/2/8.314/T)))*dv

def norm(freq, num,N_mol):
    return (N_f)*num/np.sum(num)

def F_norm(freq, num,N_mol,T,dv):
    nf=norm(freq,num,N_mol)
    return round(np.sum(8.314*T*nf*np.log(2*np.sinh(6.626*6.02*10*freq/2/8.314/T)))*dv/N_mol/1000,7)


def E_norm(freq, num,N_mol,dv):
    nf=norm(freq,num,N_mol)
    return round(np.sum(nf*dv*6.626*6.02*10*freq)/2/N_mol/1000,7)


print('degree of freedom ,6N-6:', np.sum(2*1000*out_vdos[1,1:int(out_vdos.shape[1]/2)]/8.314/10), N_f)

#a_vib = A_vib/N_mol
a_vib = F_norm(out_vdos[0,1:int(out_vdos.shape[1]/2)]/period,
       2*1000*out_vdos[1,1:int(out_vdos.shape[1]/2)]/8.314/Temp,
       N_mol,
       Temp,
       out_vdos[0,1]/period)

#pauling's residual entropy.
TSc=Temp*8.314*np.log(3/2)/1000

print('Nmol^(1/3) A_vib/Nmol(kJ/mol), TSc/Nmol/1000(kJ/mol):', np.cbrt(N_mol), a_vib, TSc)
