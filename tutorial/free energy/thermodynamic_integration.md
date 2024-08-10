# 水与冰的自由能计算

  本文以水和冰为例，分别演示了液相和固相的自由能计算方法。


## Thermodynamic Integration

### 1. Constant T (沿等温线积分)
$$
\frac{A(\rho_2,T)}{Nk_BT} = \frac{A(\rho_1,T)}{Nk_BT} + \int_{\rho_1}^{\rho_2} \frac{p(\rho)}{Nk_BT\rho^2} d\rho
$$

### 2. Constant p (沿等压线积分)
$$
\frac{G(T_2,p)}{Nk_BT_2} = \frac{G(T_1,p)}{Nk_BT_1} + \int_{T_1}^{T_2} \frac{H(T)}{Nk_BT^2} dT
$$

### 3. Constant V (沿等容线积分)
$$
\frac{A(T_2,V)}{Nk_BT_2} = \frac{A(T_1,V)}{Nk_BT_1} + \int_{T_1}^{T_2} \frac{U(T)}{Nk_BT^2} dT
$$

### 4. Hamiltonian Integration
只要在两个状态(两个状态可以指LJ流体到水 或 甲烷分子到钠离子等)间找到任意一条可逆的路径, 对该路径进行积分即可求出两个状态间的自由能变化. 此处通过定义一个耦合参数 $\lambda$ , 将两个状态的哈密顿量关联起来. 

$$
U(\lambda) = \lambda U_1 + (1-\lambda) U_2 , ~\lambda = 0\sim1
$$

$$
A(N,V,T,\lambda) = -k_BT\ln \Big{[} \frac{q^N}{N!} \int e^{-\beta U(\lambda)} d1...dN \Big{]}
$$

$$
\frac{\partial A(N,V,T,\lambda)}{\partial \lambda} =   \Big{\langle} \frac{\partial U(\lambda)}{\partial \lambda}  \Big{\rangle}_{N,V,T,\lambda}
$$

$$
A(N,V,T,\lambda=1) = A(N,V,T,\lambda=0)  +\int_{\lambda=0}^{\lambda=1}  \Big{\langle} \frac{\partial U(\lambda)}{\partial \lambda}  \Big{\rangle}_{N,V,T,\lambda} d\lambda
$$

其中q为粒子的配分函数, 以上公式在使用过程中应注意, 在选择积分路径时不应跨越相边界.

## Einstein Crystal Method

Einstein Crystal Method专门用于计算固体的自由能( $A_{sol}$ ), 与Hamiltonian Integration类似, 其将其中一个状态选定为Einstein Crystal. 该状态下, 每个粒子间无相互作用, 但每个粒子都被一个谐振势约束在晶格上, 此时晶体的自由能( $A_{Ein}$ )具有解析形式.
由此可知, 当通过解析方式求得 $A_{Ein}$ 后, 进而可由Hamiltonian Integration求得 $A_{sol}$ . 其路径如下: 

$$
\Large
\begin{align*}  
 \underbrace{ A_{Ein} \stackrel{\Delta A_1}{\longrightarrow} A_{Ein}^{CM}  \stackrel{\Delta A_2}{\longrightarrow} A_{sol}^{CM}  \stackrel{\Delta A_3}{\longrightarrow} A_{sol} } _{\Delta A= A _{sol} - A _{Ein}}
\end{align*}  
$$

其中, 与各项自由能相对应的状态为:  
$A_{Ein}^{CM}$ : 固定体系质心的Einstein Crystal  
$A_{sol}^{CM}$ : 固定体系质心的固体  
$A_{sol}$ : 真实的固体  
因此, $A_{sol}$ 可表示为:  

$$
\Large
\begin{align*}  
A_{sol} = \Delta A + A_{Ein} = \Delta A_3 +\Delta A_2 + \Delta A_1 + A_{Ein}
\end{align*}  
$$

### 1.初始状态Einstein Crystal
Einstein Crystal的Hamiltonian为:

$$
\Large
\begin{align*}  
H_{Ein} = \sum _{i} \frac{\vec{p} _i ^2}{2m _i} +\sum _{i} k(\vec{r} _i - \vec{r} _i^0)^2
\end{align*}  
$$

配分函数 $Q_{Ein}$ 为:

$$
\Large
\begin{align*}  
Q_{Ein} &= \frac{1}{N!} \underbrace{\frac{1}{h^{3N}} \int e^{-\beta \sum \frac{\vec{p} _i ^2}{2m_i} } \mathrm{d} \vec{p} _1 \cdots \mathrm{d} \vec{p} _N } _{\rm{ideal \  gas \  term}: P _{id}}  \underbrace{\int e^{-\beta \sum k(\vec{r} _i - \vec{r} _i^0)^2 } \mathrm{d} \vec{r}_1 \cdots \mathrm{d} \vec{r} _N } _{\rm{configuration \  integral}: Z}   \\
Q _{Ein} &= P _{id} \frac{1}{N!} Z
\end{align*}  
$$

其中 $\vec{r}, \vec{p}, k, h$ 分别为粒子的坐标, 动量, 谐振式的力常数和普朗克常数; $\beta=1/k_BT$ .  
对于给定的晶体结构, 每个粒子的位置并不是指定好的(例如, (x0,y0,z0)这个点可能被粒子1占据,也可能被粒子N占据, 总共有 $N!$ 种可能, 尽管晶体结构不变, 但这些结构在构型空间(configuration space)中的位置却不同), 即构型积分(configuration integral)中包含 $N!$ 种相同的排列方式, 每一种的构型积分均为 $Z_1$ , 因此 $Z= N! Z_1$ , 进而: 

$$
\Large
\begin{align*}  
Q_{Ein} =  P _{id} Z_1  
\end{align*}  
$$

若体系为单质体系, 则 $m_i = m;\ Z_1 = (\pi / \beta k )^{3N/2}$ . 自由能 $A_{Ein}$ 为:

$$
\Large
\begin{align*}  
A_{Ein} =  -\frac{\mathrm{ln} Q _{Ein}}{\beta} 
\end{align*}  
$$

### 2.由Einstein Crystal到固定质心的Einstein Crystal
#### 固定质心的必要性
可以避免数值积分的不连续性, 见文献:
New Monte Carlo method to compute the free energy of arbitrary solids. Application to the fcc and hcp phases of hard spheres

通过引入delta函数可以将固定质心后的配分函数 $Q_{Ein}^{CM}$ 表达为:  

$$
\Large
\begin{align*}  
Q_{Ein}^{CM} =  P _{id}^{CM} Z_1^{CM}  =P _{id}^{CM} \int e^{-\beta \sum k(\vec{r} _i - \vec{r} _i^0)^2 }  \delta(\underbrace{ \frac{\sum m _i (\vec{r} _i- \vec{r} _i^0)}{\sum m_i} } _{\epsilon}) \mathrm{d} \vec{r}_1 \cdots \mathrm{d} \vec{r} _N 
\end{align*}  
$$

上式中 $\epsilon = \frac{\sum m _i (r _i-r _0)}{\sum m_i}$ 为质心与固定位置之差. $P _{id}^{CM}$ 和构型积分项均未知, 因此需要从 $Q _{Ein}$ 中寻求解法.  
对于 $Q _{Ein}$ , 相应地, 可以做同样的等价变换:

$$
\Large
\begin{align*}  
Q_{Ein} &=  P _{id} Z_1  =P _{id} \int e^{-\beta \sum k(\vec{r} _i - \vec{r} _i^0)^2 }  \delta( \underbrace {\epsilon - \frac{\sum m _i (\vec{r} _i-\vec{r} _i^0)}{\sum m_i} } _{\equiv 0}) \mathrm{d} \vec{r}_1 \cdots \mathrm{d} \vec{r} _N \mathrm{d} \epsilon \\
&= P _{id} \int e^{-\beta \sum k(\hat{r} _i - \vec{r} _i^0 + \epsilon)^2 }  \delta( \underbrace { \frac{\sum m _i (\hat{r} _i- \vec{r} _i^0)}{\sum m_i} } _{changing\  variables: \hat{r} _i = \vec{r} _i -\epsilon}) \mathrm{d} \hat{r}_1 \cdots \mathrm{d} \hat{r} _N \mathrm{d} \epsilon \\
&= P _{id} \int e^{-\beta \sum k(\hat{r} _i - \vec{r} _i^0)^2}  e^{-\beta N k\epsilon ^2 }  \delta( \frac{\sum m _i (\hat{r} _i- \vec{r} _i^0)}{\sum m_i}) \mathrm{d} \hat{r}_1 \cdots \mathrm{d} \hat{r} _N \mathrm{d} \epsilon \\
&= P _{id} \int e^{-\beta N k\epsilon ^2 } \mathrm{d} \epsilon \int e^{-\beta \sum k(\hat{r} _i - \vec{r} _i^0)^2}    \delta( \frac{\sum m _i (\hat{r} _i- \vec{r} _i^0)}{\sum m_i}) \mathrm{d} \hat{r}_1 \cdots \mathrm{d} \hat{r} _N  \\
&= P _{id} (\frac{2 \pi}{k \beta})^{3/2}  \int e^{-\beta \sum k(\hat{r} _i - \vec{r} _i^0)^2}    \delta( \frac{\sum m _i (\hat{r} _i- \vec{r} _i^0)}{\sum m_i}) \mathrm{d} \hat{r}_1 \cdots \mathrm{d} \hat{r} _N  \\
\end{align*}  
$$

对比 $Q_{Ein}$ 与 $Q_{Ein}^{CM}$ 可以发现积分项完全相同, 因此可以消除积分项:

$$
\Large
\begin{align*}  
\frac{Q_{Ein} }{Q_{Ein}^{CM}} = \frac{P _{id}}{P _{id} ^{CM}} (\frac{2 \pi}{k \beta})^{3/2}
\end{align*}  
$$

自然地, 可知自由能差 $\Delta A_1$ 为:

$$
\Large
\begin{align*}  
-\Delta A_1 &= A _{Ein} - A _{Ein}^{CM} = -\frac{1}{\beta}\mathrm{ln}(\frac{P _{id}}{P _{id} ^{CM}}) -\frac{1}{\beta} \mathrm{ln} (\frac{2 \pi}{k \beta})^{3/2} \\
\Delta A_1 &= \frac{1}{\beta}\mathrm{ln}(\frac{P _{id}}{P _{id} ^{CM}}) +\frac{1}{\beta} \mathrm{ln} (\frac{2 \pi}{k \beta})^{3/2}
\end{align*}  
$$

### 3.由固定质心的Einstein Crystal到固定质心的真实固体
在固定质心的状态下 $H_{Ein}$ 和 $H_{sol}$ 分别为:

$$
\Large
\begin{align*}  
H _{Ein} &= \sum _{i} \frac{\vec{p} _i ^2}{2m _i} +\sum _{i} k(\vec{r} _i - \vec{r} _i^0)^2 \\
H _{sol} &= \sum _{i} \frac{\vec{p} _i ^2}{2m _i} +\sum U _{sol}(r _{ij})
\end{align*}  
$$

对两个状态做Hamiltonian integration即可求得 $\Delta A_2$ .

### 4.由固定质心的真实固体到真实固体
与第二步相同, 从配分函数入手, 如法炮制: 

$$
\Large
\begin{align*}  
Q _{sol} &= \frac{P _{id}}{N!} \int e^{-\beta U _{sol}(r) }  \mathrm{d} \vec{r}_1 \cdots \mathrm{d} \vec{r} _N \\
Q _{sol}^{CM} &= \frac{P _{id}^{CM}}{N!} \int e^{-\beta U _{sol}(r) }  \delta( \underbrace {\frac{\sum m _i (\vec{r} _i)}{\sum m_i}} _{COM\ at\ the\ origin})  \mathrm{d} \vec{r}_1 \cdots \mathrm{d} \vec{r} _N 
\end{align*}  
$$

真实固体与Einstein Crystal的区别是前者的势能由原子间的相对距离决定,而后者则由原子相对于其平衡位置的距离决定; 由此可知, 任意平移真实固体, 其势能均不会发生改变(即 $U(r) = U(r-r0)$ ). 将delta函数引入真实固体的配分函数, 假设对固体平移到任意位置 $r_0$ :

$$
\Large
\begin{align*}  
Q _{sol} &= \frac{P _{id}}{N!} \int e^{-\beta U _{sol}(r-r0) }  \delta( r_0 -\frac{\sum m _i (\vec{r} _i)}{\sum m_i} )  \mathrm{d}r _0 \mathrm{d} \vec{r}_1 \cdots \mathrm{d} \vec{r} _N \\
&= \frac{P _{id}}{N!} \int e^{-\beta U _{sol}(r) }  \delta(\frac{\sum m _i (\vec{r} _i-r_0)}{\sum m_i} )  \mathrm{d}r _0 \mathrm{d} \vec{r}_1 \cdots \mathrm{d} \vec{r} _N \\
&= \frac{P _{id}}{N!} \int \mathrm{d}r _0 \int e^{-\beta U _{sol}(r) }  \delta(\frac{\sum m _i (\vec{r})}{\sum m_i} )  \mathrm{d} \vec{r}_1 \cdots \mathrm{d} \vec{r} _N \\
&= \frac{P _{id}}{N!} V \int e^{-\beta U _{sol}(r) } \delta(\frac{\sum m _i (\vec{r})}{\sum m_i} ) \mathrm{d} \vec{r}_1 \cdots \mathrm{d} \vec{r} _N \\
&= \frac{P _{id}}{N!} V (N-1)! \int_1 e^{-\beta U _{sol}(r) } \delta(\frac{\sum m _i (\vec{r})}{\sum m_i} ) \mathrm{d} \vec{r}_1 \cdots \mathrm{d} \vec{r} _N \\
&= \frac{P _{id}V}{N} \int_1 e^{-\beta U _{sol}(r) } \delta(\frac{\sum m _i (\vec{r})}{\sum m_i} ) \mathrm{d} \vec{r}_1 \cdots \mathrm{d} \vec{r} _N \\
\end{align*}  
$$

在上式中, 积分项里最终固定了质心, 因此自由度减少为了N-1, 因此其只具有 (N-1)! 种等同的排列方式.  对比 $Q_{sol}$ 与 $Q_{sol}^{CM}$ 可以发现积分项完全相同, 因此可以消除积分项:

$$
\Large
\begin{align*}  
\frac{Q_{sol} }{Q_{sol}^{CM}} = \frac{P _{id}}{P _{id} ^{CM}} (\frac{V}{N})
\end{align*}  
$$

自然地, 可知自由能差 $\Delta A_3$ 为:

$$
\Large
\begin{align*}  
\Delta A_3 &= A _{sol} - A _{sol}^{CM} = -\frac{1}{\beta}\mathrm{ln}(\frac{P _{id}}{P _{id} ^{CM}}) +\frac{1}{\beta} \mathrm{ln} \rho \\
\end{align*}  
$$

可以注意到, $\frac{1}{\beta}\mathrm{ln}(\frac{P _{id}}{P _{id} ^{CM}})$ 可以与 $\Delta A_1$ 中的项相消, 进而避免了求解未知项 $P _{id} ^{CM}$ 的问题.


### Einstein Molecule Method
Einstein Crystal方法中, 为了避免数值积分不连续的问题, 采用了固定质心的方法. 除此之外, 这一问题也可以通过固定一个粒子的位置(例如固定粒子1的位置, 粒子1依然可以自由转动)来解决, 即Einstein Molecule法. 
类比Einstein Crystal可以写出其Hamiltonian, 由于粒子1的位置被固定在了平衡位置, 因此其势能项从2开始计算.

$$
\Large
\begin{align*}  
H _{Ein-mol} &= \sum _{i=1} \frac{\vec{p} _i ^2}{2m _i} +\sum _{i=2} k(\vec{r} _i - \vec{r} _i^0)^2 
\end{align*}  
$$



## 水的自由能计算[^1]

  液态水自由能的绝对值可以通过Hamiltonian integration进行计算. 其中始末状态分别选为LJ流体以及水即可. 目前主流的水模型中均只有氧原子有LJ参数, 氢原子为点电荷. 因此只需要逐渐关闭水分子的静电作用, 即减小原子电荷, 便可实现从LJ流体到水的可逆变化. 该过程的自由能如下:

$$
U(\lambda) = \lambda U_{TIP4P} + (1-\lambda) U_{LJ} , ~\lambda = 0\sim1
$$

$$
 = \lambda (U_{Coul} + U_{LJ}) + (1-\lambda) U_{LJ} 
$$
<!--
$$
 = \lambda (k\frac{q_1 q_2}{r^2} + U_{LJ}) + (1-\lambda) U_{LJ} 
$$

$$
 =  k\frac{\sqrt{\lambda}q_1 \sqrt{\lambda}q_2}{r^2} + \lambda(U_{LJ}) + (1-\lambda) U_{LJ} 
$$
-->
$$
A_{water}(N,V,T) - A_{LJ}(N,V,T)  = \int_{\lambda=0}^{\lambda=1}  \Big{\langle} \frac{\partial U(\lambda)}{\partial \lambda} \Big{\rangle}_{N,V,T,\lambda} d\lambda 
$$

$$
=\int_{\lambda=0}^{\lambda=1}  \Big{\langle} U_{water} - U_{LJ} \Big{\rangle}_{N,V,T,\lambda} d\lambda
$$

$$
=\int_{\lambda=0}^{\lambda=1}  \Big{\langle} U_{Coul} \Big{\rangle}_{N,V,T,\lambda} d\lambda
$$

其中, $U_{Coul}$ 为体系的静电势能, 任意温度压力下, LJ流体的自由能可以通过状态方程算出(见MBWR.py[^2]). 以上公式的含义为: 首先在不同的 $\lambda$ 下进行MD模拟产生模拟轨迹, 随后对这些轨迹计算 $\lambda=1$ 下的静电势能, 最后求积分即可. 在MD模拟中, 逐步关闭静电作用是通过减小原子电荷实现的, 即给电荷乘一个缩放因子, 由库仑定律可知, 电荷的缩放因子应为 $\sqrt{\lambda}$.

$$
\lambda U_{Coul} = k\frac{\sqrt{\lambda}q_i \sqrt{\lambda}q_j}{r^2}
$$

下文中将采用SPC/E水模型为例演示自由能计算过程, 相关文件均可在water路径下获得. 
SPC/E水模型的参数如下:

|  $\epsilon /k_B$ (K)   | $\sigma$ (nm)  | $q_H$ (e) |
|  ----                | ----         | ----    |
| 78.20                | 0.31656      | 0.4238  |

#### 1. Equilibration
   
   创建一个360个水分子的水盒子, 在225K 564bar下进行1ns的MD模拟, 将体系跑平衡. 
   
#### 2. Hamitonian Integration
   
   模拟流程如下: 在11个 $\lambda$ 下分别进行能量最小化, 平衡相以及产生相模拟. 随后使用SPC/E的原子电荷对产生相中的11段轨迹进行rerun, 计算静电势能, 最终积分求解自由能. 
   
   初始结构: 以步骤1中产生的结构作为初始结构 (注意:此处与[^1]中的操作不同, 文献中采用前一 $\lambda$ 下产生的结构作为后一 $\lambda$ 下的初始结构. 因为平衡态下的模拟结果不应依赖于初始结构, 所以这一差异对结果没有影响.)
 
   拓扑文件: 以0.1为间隔从0到1产生11个 $\lambda$, 对水分子的原子电荷进行放缩( $q_H(\lambda) = \sqrt{\lambda} q_H$ ), 并写入拓扑文件中.

   模拟参数: 模拟均在NVT系综下进行, 非键作用的cut-off和rlist均设为0.85, 长程静电作用通过PME方法计算, 水分子的键长键角均要进行约束.

   文献值如下[^1] (此处德布罗意热力学波长 $\Lambda$ = 0.1 nm, $A_{LJ} =A_{LJ}^{id} + A_{LJ}^{res}$):
   
| Model   | $T$ (K)  | $p$ (bar) | $\rho$  (g/cm $^3$ ) | $A_{water}/(Nk_BT)$ | $A_{LJ}^{res}/(Nk_BT)$ | $A_{LJ}^{id}/(Nk_BT)$ | $\Delta A/(Nk_BT)$ |
|  ----   | ----     | ----      | ----              | ----              | ----                   | ----                  | ----               | 
| SPC/E | 225 | 564 | 1.05 | -21.82 | 2.500 | -4.350 | -19.97|

其中 $A_{LJ}^{id} = Nk_BT(\ln(\rho\Lambda^3)-1)$

<!--
```python
    #!/usr/bin/env python3
    print("Hello, World!");
```
-->


## 冰的自由能计算[^2]


## 参考文献
[^1]: J. Phys.: Condens. Matter 20 (2008) 153101
[^2]:  Johnson J K, Zollweg J A and Gubbins K E 1993 Mol. Phys. 78-591. The Lennard-Jones equation of state revisited.

