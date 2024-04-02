# Mold Integration for interfacial free energy[^1]

在熔点下, 水中形成一层冰晶面所产生的自由能变化为 $\Delta G = 2A\gamma$,  $A$ 为冰晶面的面积, 由于晶面上下两侧均与水接触, 因此冰水界面面积为 $2A$. 
由此可求出界面自由能为 $\gamma = \Delta G/ 2A$ . 由于均相成核十分困难, 因此需要使用一个势阱诱导冰晶面的形成; 为了使水分子按晶格排列, 这个势阱需要足够窄; 通过逐渐开启势阱与水分子的相互作用, 晶面的形成自由能可以通过热力学积分获得.

## 热力学积分:
由于存在势阱, 因此热力学积分求得的晶面形成自由能 $\Delta G^m$ 不仅是水中形成晶面过程的自由能变化 $\Delta G$ , 实际上还包含了势阱部分的贡献 ( $-N_w \epsilon$ )

体系势能 $U(\lambda)$ 为: 

$$
U(\lambda)= U_{pp} ( \mathbf{r_1, ..., r_N} )  + \lambda U_{pm}( \mathbf{r_1, ..., r_N ; r_{w1}, ..., r_{wN_w}})
$$

$$
U_{pm} = \sum_{i=1}^{i=N} \sum_{w_j=1}^{w_j=N_w} u_{pw}(r_{iw_j})
$$

其中 $N$ -原子数, $N_w$ -势阱个数, $U_{pp}$ -原子间势能, $U_{pm}$ -势阱和原子间的势能, $\lambda = 0 \sim 1$ -耦合常数, $u_{pw}$ 是势阱和粒子间的相互作用(矩形势: $u_{pw} = -\epsilon , r_{iw_j}<=r_w ; =0 , otherwise$ , $r_w$ 和 $\epsilon$ 为矩形势的宽度和深度.)

通过热力学积分可求得:

$$
\Delta G^m = \int_{\lambda=0}^{\lambda=1}  \big{\langle} \frac{\partial U(\lambda)}{\partial \lambda}  \big{\rangle} _{\lambda, N, p_z, T} d\lambda
$$

$$
\Delta G^m = \int_{\lambda=0}^{\lambda=1}  \big{\langle}  U_{pm}  \big{\rangle} _{\lambda, N, p_z, T} d\lambda
$$

$p_z$ 表示只对垂直于晶面的方向控压


[^1]:J. Chem. Phys. 141, 134709 (2014)
[^2]:J. Phys. Chem. C 2016, 120, 8068−8075
