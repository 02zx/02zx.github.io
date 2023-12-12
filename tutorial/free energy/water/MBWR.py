import numpy as np

class reduced_unit:
    #Usage: reduced_unit(300, 360, 2.71**3).cal_T() 
    def __init__(self,T,N,V,epsilon=78.20*8.314/1000, sigma=0.31656):
        #The default parameters are for SPC/E water model.
        self.T=T
        self.N=N
        self.V=V
        self.epsilon=epsilon
        self.sigma=sigma
        
    def h(self):
        print("T-K, V-nm^3, epsilon-kJ/mol, sigma-nm")
    
    def cal_T(self):
        return 8.314*self.T/(self.epsilon*1000)
    
    def cal_rho(self):
        return (self.N*self.sigma**3)/self.V
    
    

def MBWR(T,rho):
    #T and rho are in reduced units.

    #Johnson J K, Zollweg J A and Gubbins K E 1993 Mol.Phys. 78 591
    #See Tab.3 for validation
    
    #Tab.10
    x1=0.8623085097507421
    x2=2.976218765822098
    x3=-8.402230115796038
    x4=0.1054136629203555
    x5=-0.8564583828174598
    x6=1.582759470107601
    x7=0.7639421948305453
    x8=1.753173414312048
    x9=2.798291772190376e+3
    x10=-4.8394220260857657e-2
    x11=0.9963265197721935
    x12=-3.698000291272493e+1
    x13=2.084012299434647e+1
    x14=8.305402124717285e+1
    x15=-9.574799715203068e+2
    x16=-1.477746229234994e+2
    x17=6.398607852471505e+1
    x18=1.603993673294834e+1
    x19=6.805916615864377e+1
    x20=-2.791293578795945e+3
    x21=-6.245128304568454
    x22=-8.116836104958410e+3
    x23=1.488735559561229e+1
    x24=-1.059346754655084e+4
    x25=-1.131607632802822e+2
    x26=-8.867771540418822e+3
    x27=-3.986982844450543e+1
    x28=-4.689270299917261e+3
    x29=2.593535277438717e+2
    x30=-2.694523589434903e+3
    x31=-7.218487631550215e+2
    x32=1.721802063863269e+2
    
    #Tab. 5
    a1=x1*T+x2*np.sqrt(T)+x3+x4/T+x5/T**2
    a2=x6*T+x7+x8/T+x9/T**2
    a3=x10*T+x11+x12/T
    a4=x13
    a5=x14/T+x15/T**2
    a6=x16/T
    a7=x17/T+x18/T**2
    a8=x19/T**2
    
    #Tab.6
    b1=x20/T**2+x21/T**3
    b2=x22/T**2+x23/T**4
    b3=x24/T**2+x25/T**3
    b4=x26/T**2+x27/T**4
    b5=x28/T**2+x29/T**3
    b6=x30/T**2+x31/T**3+x32/T**4
    
    #Tab.7
    gamma=3
    F=np.exp(-gamma*rho**2)
    G1=(1-F)/(2*gamma)
    G2=-(F*rho**2-2*G1)/(2*gamma)
    G3=-(F*rho**4-4*G2)/(2*gamma)
    G4=-(F*rho**6-6*G3)/(2*gamma)
    G5=-(F*rho**8-8*G4)/(2*gamma)
    G6=-(F*rho**10-10*G5)/(2*gamma)
    
    #Tab.8
    c1=x2*np.sqrt(T)/2+x3+2*x4/T+3*x5/T**2
    c2=x7+2*x8/T+3*x9/T**2
    c3=x11+2*x12/T
    c4=x13
    c5=2*x14/T+3*x15/T**2
    c6=2*x16/T
    c7=2*x17/T+3*x18/T**2
    c8=3*x19/T**2
    
    #Tab.9
    d1=3*x20/T**2+4*x21/T**3
    d2=3*x22/T**2+5*x23/T**4
    d3=3*x24/T**2+4*x25/T**3
    d4=3*x26/T**2+5*x27/T**4
    d5=3*x28/T**2+4*x29/T**3
    d6=3*x30/T**2+4*x31/T**3+5*x32/T**4
    


        

    #Ar(NVT)=A(NVT)-Aid(NVT); Ar_ = Ar/N/epsilon
    #eq.5
    Ar_=(a1*rho**1)/1 \
      + (a2*rho**2)/2 \
      + (a3*rho**3)/3 \
      + (a4*rho**4)/4 \
      + (a5*rho**5)/5 \
      + (a6*rho**6)/6 \
      + (a7*rho**7)/7 \
      + (a8*rho**8)/8 \
      + b1*G1 + b2*G2 + b3*G3 + b4*G4 + b5*G5 + b6*G6
    

    #eq.7 P_ = P*sigma**3/epsilon
    P_=rho*T \
     + a1*rho**2 + a2*rho**3 + a3*rho**4 + a4*rho**5 + a5*rho**6 + a6*rho**7 + a7*rho**8 \
     + F*(b1*rho**3 + b2*rho**5 + b3*rho**7 + b4*rho**9 + b5*rho**11)
    
    
    #Ur_ = Ur/N/epsilon
    #eq.9 
    Ur_=(c1*rho**1)/1 \
      + (c2*rho**2)/2 \
      + (c3*rho**3)/3 \
      + (c4*rho**4)/4 \
      + (c5*rho**5)/5 \
      + (c6*rho**6)/6 \
      + (c7*rho**7)/7 \
      + (c8*rho**8)/8 \
      + d1*G1 + d2*G2 + d3*G3 + d4*G4 + d5*G5 + d6*G6

    
    #eq.10
    Gr_ = Ar_ + P_/rho - T

    return Ar_, P_, Ur_, Gr_

print('A_res, P_reduced, U_res, G_res:',
      MBWR(reduced_unit(300, 360, 2.71**3).cal_T(),reduced_unit(300, 360, 2.71**3).cal_rho()))
    
