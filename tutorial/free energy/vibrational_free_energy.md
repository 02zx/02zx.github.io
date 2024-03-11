# 计算冰在低温下的自由能[^1]

## 理论部分

当温度足够低时, 晶体中的原子在平衡位置附近振动, 这时系统的势能 $E$ 可以在平衡位置处展开成原子位移的幂级数. 

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

由此可写出正则系综下的配分函数 $Q(\beta, N, V)$ ,此处取 $U_0=0$ :

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

如果已知声子态随频率的分布 $g(\nu)$ ，又称为态密度, 则可以从配分函数求出自由能 $A$ :

$$
A = -\frac{1}{\beta} \ln Q = \frac{1}{\beta} \int_{0}^{\infty} d \nu g(\nu) \ln \bigg{\[} 2\sinh (\frac{\beta h\nu_l}{2})  \bigg{\]}
$$

## 正则分析

通过正则分析可以获得态密度, 目前gromacs对三点水模型的正则分析支持良好, 但四点水模型中正则分析会出现大量虚频无法解决。

## 速度自相关函数求算态密度

除正则分析外还可以通过对质量权重的速度自相关函数(mass-weighted velocity autocorrelation function, 简写为mw-vacf)进行傅里叶变换来获得速度的态密度分布(VDOS), 在低温下其结果与正则分析基本吻合. 通过[^2]中的公式22-23可以计算体系的自由度, 其数量应为3N, 其中N为体系中的总原子数. 实际计算出的结果会与3N略有差异.

需要注意以下几点:
1. gmx velacc只能计算vacf, 若要使用公式22需要对每种原子进行计算然后乘以其质量后求和. 另外, 按此步骤算出的mw-vacf与gmx dos的结果不同, 暂时无法确定gmx dos是如何计算mw-vacf的. 同时需要注意, gmx dos计算的VDOS可能存在问题[^3]
2. 该方法基于谐振近似, 应在尽可能低的温度下进行模拟
3. 体积对计算结果又显著影响,若要和正则分析的结果对比应保持相同体积(NVT下进行模拟).



[^1]:Introduction to Modern Statistical Mechanics-pp:90-92
[^2]:J. Chem. Phys. 150, 194111 (2019)
[^3]:Two Faces of the Two-Phase Thermodynamic Model. J.Chem.TheoryComput.2021, 17, 7187−7194
[^4]:J. Phys. Chem. B2010,114,8191–8198
[^5]:Nanoscale, 2020, 12, 18701
