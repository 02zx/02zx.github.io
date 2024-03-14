# 计算冰在低温下的自由能[^1]

## 理论部分

当温度足够低时, 晶体中的 $N$ 个原子在平衡位置附近振动, 这时系统的势能 $E$ 可以在平衡位置处展开成原子位移的幂级数. 

$$
E=U_0 + \sum_{i=1}^{N}\sum_{a=x,y,z} (S_{i,a} - S_{i,a}^0) K_i,a + \frac{1}{2} \sum_{i,j=1}^{N}\sum_{\alpha, \gamma=x,y,z} (S_{i,\alpha}-S_{i,\alpha}^0) (S_{j,\gamma}-S_{j,\gamma}^0) K_{i,\alpha,j,\gamma} + ...
$$

其中 $U_0$ 是零点能(Zero point energy)-势能的极小点, $S_{i,a}$ 是粒子 $i$ 的 $a$ 轴坐标,  $S_{i,a}^0$ 是相应的平衡位置. 因此二次项系数 $K_{i,\alpha,j,\gamma}$ 具有力常数的单位.
此外, 式中一次项总是为0, 因为极小点处的一阶导总是为0.

在二阶项处截断的情况称之为谐振近似:

$$
E=U_0 + \frac{1}{2} \sum_{i,j=1}^{N}\sum_{\alpha, \gamma=x,y,z} (S_{i,\alpha}-S_{i,\alpha}^0) (S_{j,\gamma}-S_{j,\gamma}^0) K_{i,\alpha,j,\gamma} 
$$

上式中的 $K_{i,\alpha,j,\gamma}$ 是一个DN $\times$ DN 阶矩阵, D为坐标维度(通常是3维空间笛卡尔坐标,x,y,z).

$$
 K_{i,\alpha,j,\gamma} = \begin{bmatrix}  
    K_{1,x,1,x} & K_{1,x,1,y} &  K_{1,x,1,z} & K_{1,x,2,x} & ... & K_{1,x,n,z}\\  
    K_{1,y,1,x} & K_{1,y,1,y} &  K_{1,y,1,z} & K_{1,y,2,x} & ... & K_{1,y,n,z}\\  
    K_{1,z,1,x} & K_{1,z,1,y} &  K_{1,z,1,z} & K_{1,z,2,x} & ... & K_{1,z,n,z}\\  
    K_{2,x,1,x} & K_{2,x,1,y} &  K_{2,x,1,z} & K_{2,x,2,x} & ... & K_{2,x,n,z}\\  
    ...         & ...         & ...          &...          & ... & ...         \\  
    K_{n,z,1,x} & K_{n,z,1,y} &  K_{n,z,1,z} & K_{n,z,2,x} & ... & K_{n,z,2,z}  \\
\end{bmatrix}  
$$

以上矩阵是对称的, 因此可解出DN个一维的normal modes, 每个normal mode都是独立的. 所以势能(其哈密顿量为 $H$ )也就可以表示成DN个normal modes组成的二次函数之和( $H_{l}$ ):

$$
H=\sum_{l=1}^{DN} H_{l} 
$$

谐振子的本征能量为 $(\frac{1}{2} + n)h\nu, n=0,1,2,...$ . (频率 $\nu$ 和 力常数 $K$ 的关系是 $\nu=\frac{1}{2\pi}\sqrt{K/\mu}$, $\mu$ 为约化质量 )
因此上式中每个normal modes都对应着一个基频 $\nu_l$ . 

此处, 将振动频率为 $\nu_l$ 的谐振子称为声子, $n_l$ 为处于能量 $h\nu_l$ 的声子态的声子数.


$$
E= \sum_{l=1}^{DN} \frac{1}{2} h\nu_l + \sum_{l=1}^{DN} n_l h\nu_l + U_0
$$

由此可写出正则系综下的配分函数 $Q(\beta, N, V)$ ,此处取 $U_0=0$, 这时自由能全部来自于谐振项, 此处将其记为 $A_{vib}$ :

$$
Q(\beta, N, V) = \sum_{n_1,n_2,...,=0}^{\infty} e^{-\beta \sum_{l} (\frac{1}{2}+n_l) h\nu_l}
$$

$$
Q(\beta, N, V) = \prod_{l=1}^{DN}  \sum_{n} e^{-\beta (\frac{1}{2}+n) h\nu_l}
$$

$$
\ln(Q) = -\sum_{l=1}^{DN} ln(e^{\beta h\nu_l/2} - e^{-\beta h \nu_l/2}) ; \sum_k^{\infty} a r^k = \frac{a}{1-r}, |k|<1
$$

$$
\ln(Q) = -\sum_{l=1}^{DN} \ln \bigg{\[} 2\sinh (\frac{\beta h\nu_l}{2})  \bigg{\]}
$$

如果已知声子态随频率的分布 $g(\nu)$ ，又称为态密度, 则可以从配分函数求出自由能 $A_{vib}$ :

$$
A_{vib} = -\frac{1}{\beta} \ln Q = \frac{1}{\beta} \int_{0}^{\infty} d \nu g(\nu) \ln \bigg{\[} 2\sinh (\frac{\beta h\nu_l}{2})  \bigg{\]} ~~~~~~~Eq.1
$$

总自由能 $A$ 为:

$$
A = U_0(V) + A_{vib}(T,V) -TS_c ~~~~~~~Eq.2
$$

其中 $S_c = k_B N \ln (3/2)$ 称为 Pauling's residual entropy.


## 正则分析

通过正则分析可以获得态密度, 目前gromacs对三点水模型的正则分析支持良好, 但四点水模型中正则分析会出现大量虚频无法解决。

## 速度自相关函数求算态密度

除正则分析外还可以通过对质量权重的速度自相关函数(mass-weighted velocity autocorrelation function, 简写为mw-vacf)进行傅里叶变换来获得速度的态密度分布(VDOS)[^2], 在低温下其结果与正则分析基本吻合. 

$$
g(\nu) = 2\beta \sum_{j=1}^N m_j \int dt \langle \vec{v}_j (\tau) \vec{v}_j (\tau + t) \rangle _\tau e^{-i2\pi \nu t} ~~~~~~~Eq.3
$$

通过对态密度进行积分可以获得体系的自由度( $N_f$ ), 其数量应为3N (对于含有约束的体系应另作考虑),  其中N为体系中的总原子数. $N_f$ 可用于验证结果的准确性, 通常实际计算出的结果会与3N略有差异.

$$
N_f = \int_0^{\infty} d\nu g(\nu) ~~~~~~~Eq.4
$$



## 算例
水模型:TIP4P/ICE-rigid, 模拟体系: 冰Ih
### 计算速度自相关函数
#### 体系预平衡
准备好Ih的gro文件后首先对其进行能量最小化, 随后在10K下进行 $NpT$ 模拟(若是模拟团簇则用 $NVT$ ), mdp中注意事项如下:

```bash
#---npt.mdp----
constraints         =  hbonds
dt                  =  0.001
tcoupl              =  V-rescale
tau-t               =  0.1
ref-t               =  10
Pcoupl              =  C-rescale
Pcoupltype          =  isotropic
tau_p               =  1
compressibility     =  4.5e-5
ref_p               =  1.0
#以下用于团簇模拟, 晶体不需要设置
comm-mode           = Angular
```
#### $NVE$ 模拟
模拟时间和步长可能会影响最终结果, 需要进行测试

在 $NVE$ 系综下模拟50ps用于计算自相关函数
```bash
#---nve.mdp----
constraints         =  hbonds
dt                  =  0.001
nsteps              =  50000
nstvout             =  1
continuation        = yes
gen_vel             = no
#以下用于团簇模拟, 晶体不需要设置
comm-mode           = Angular
```
分析自相关函数前首先要将不同类型的原子分类做成index.ndx文件

```bash
mkdir ndx
gmx_mpi make_ndx -f nve/$i.gro -o ./ndx/${i}.ndx << EOF
 a OW
 a HW1 HW2
q
EOF
```

随后使用gmx velacc对每种类型的原子进行分析

```bash
mkdir $i
echo "3"|gmx_mpi velacc -f ./md/$i.trr -s ./md/$i.tpr  -o ./$i/O.xvg -n ./ndx/${i}.ndx -dt 0.01 -b 10 -acflen 501 -nonormalize -norecip
echo "4"|gmx_mpi velacc -f ./md/$i.trr -s ./md/$i.tpr  -o ./$i/HW.xvg -n ./ndx/${i}.ndx -dt 0.01 -b 10 -acflen 501 -nonormalize -norecip

sed -i "1,17"d $i/O.xvg
sed -i "$"d $i/O.xvg
sed -i "1,17"d $i/HW.xvg
sed -i "$"d $i/HW.xvg
```


需要注意以下几点:
1. gmx velacc只能计算vacf, 若要使用[^2]的公式22, 需要对每种原子进行计算然后乘以其质量后求和. 另外, 按此步骤算出的mw-vacf与gmx dos的结果不同, 暂时无法确定gmx dos是如何计算mw-vacf的, gmx dos计算的VDOS可能存在问题[^3]
2. 该方法基于谐振近似, 应在尽可能低的温度下进行模拟
3. 体积对计算结果又显著影响,若要和正则分析的结果对比应保持相同体积(NVT下进行模拟).

### 计算态密度 $g(\nu)$ 
计算mw-vacf:

```python
#import matplotlib.pyplot as plt
import numpy as np
import sys

N_mol=int(sys.argv[1])  # number of molecules sys.argv[1]
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
```

对mw-vacf进行傅里叶变换:

```python
def vdos(N_mol):
    filename = "./{}/water.mwvacf".format(N_mol)
    data = np.loadtxt(filename)
    
    x = data[:, 0]
    y = data[:, 1]
    
    
    N=len(y)
    print(N)
    dt = data[1, 0]
    L = dt*float(N-1)
    
    print(L)
    
    for i in range(len(y)):
        if i == 0:
            continue
        y[i] = y[i] * np.sin(np.pi * i / (N-1)) / (np.pi * i / (N-1)) 
    
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
    print('#frequency(cm^-1) normalized_vdos, n, vdos',file=f)
    for i in range(int(out_vdos.shape[1]/2)):
        print(round(out_vdos[0,i]/period/(3/100),4),
              round(out_vdos[1,i]/total,4),
              round(out_vdos[0,i],4),
              round(out_vdos[1,i],4)
              ,file=f)

```

### 对态密度积分求得振动自由能( $A_{vib}$ )
对于rigid water model, 首先需要考虑约束, 每个水分子存在两个 $O-H$ 键的键长约束和一个 $H-H$ 的长度约束(相当于约束了HOH的角度), 因此对于 $N$ 个原子的体系共存在 $N$ 个约束;另外还需要排除整个体系在3个维度上的平动和转动自由度, 因此体系的自由度应为 $N_f = 3N-N-6= 2N -6 = 6N_{mol}-6$.



```python
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
```
### 计算总自由能( $A$ )
上一步中已经求得 $A_{vib}, TS_c$ 了, 只需再加上 $U_0$ . 其数值可以通过能量最小化获得, 使用最速下降法计算即可 (gromacs中只有该方法支持约束); 最终带入 Eq.2 即可.

## 待补充内容
Pauling's residual entropy, 1PT+AC, 2PT method

[^1]:Introduction to Modern Statistical Mechanics-pp:90-92
[^2]:J. Chem. Phys. 150, 194111 (2019)
[^3]:Two Faces of the Two-Phase Thermodynamic Model. J.Chem.TheoryComput.2021, 17, 7187−7194
[^4]:J. Phys. Chem. B2010,114,8191–8198
[^5]:Nanoscale, 2020, 12, 18701
