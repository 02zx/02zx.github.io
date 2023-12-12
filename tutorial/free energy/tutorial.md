# 水与冰的自由能计算

  本文以水和冰为例，分别演示了液相和固相的自由能计算方法。


## Thermodynamic integration

### 1. constant T (沿等温线积分)
$$
\frac{A(\rho_2,T)}{Nk_BT} = \frac{A(\rho_1,T)}{Nk_BT} + \int_{\rho_1}^{\rho_2} \frac{p(\rho)}{Nk_BT\rho^2} d\rho
$$

### 2. constant p (沿等压线积分)
$$
\frac{G(T_2,p)}{Nk_BT_2} = \frac{G(T_1,p)}{Nk_BT_1} + \int_{T_1}^{T_2} \frac{H(T)}{Nk_BT^2} dT
$$

### 3. constant V (沿等容线积分)
$$
\frac{A(T_2,V)}{Nk_BT_2} = \frac{A(T_1,V)}{Nk_BT_1} + \int_{T_1}^{T_2} \frac{U(T)}{Nk_BT^2} dT
$$

### 4. Hamiltonian integration
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

以上公式在使用过程中应注意, 在选择积分路径时不应跨越相边界.

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

下文中将采用SPC/E水模型为例演示自由能计算过程.
```python
    #!/usr/bin/env python3
    print("Hello, World!");
```



## 参考文献
[^1]: J. Phys.: Condens. Matter 20 (2008) 153101
[^2]:  Johnson J K, Zollweg J A and Gubbins K E 1993 Mol. Phys. 78-591.

