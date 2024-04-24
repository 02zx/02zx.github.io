
# The Clausius virial theorem

$$
\Large
\begin{align*}  
  m_i \ddot{r}  & = F_i \\
  \frac{\mathrm{d}}{\mathrm{d}t	} (r\dot{r})  & = \dot{r}^2 + r\ddot{r} \\
   r_i \ddot{r}_i  & = \frac{\mathrm{d}}{\mathrm{d}t	} (r_i \dot{r}_i) - \dot{r}_i^2 \\
   \overbrace{m_i r_i \ddot{r}_i}^{F_i r_i} & = \frac{\mathrm{d}}{\mathrm{d}t	} (m_i r_i \dot{r}_i)  - m_i \dot{r}_i^2 \\
 \lim _{\tau \to \infty} \langle \sum^{N} _{i} F_i r_i \rangle _\tau & =  \underbrace{ \lim _{\tau \to \infty}  \langle \sum^{N} _{i} \frac{\mathrm{d}}{\mathrm{d}t	} (m_i r_i \dot{r}_i) \rangle _\tau  } -  \lim _{\tau \to \infty}  \langle \sum^{N} _{i} m_i  \dot{r}_i^2 \rangle _\tau \\
  & {\lim _{\tau \to \infty} \sum^{N} _{i} \frac{m_i r_i(\tau) \dot{r}_i(\tau) -  m_i r_i(0) \dot{r}_i(0)}{\tau}} = 0 \\
  \underbrace{\lim _{\tau \to \infty} \langle \sum^{N} _{i} F_i r_i \rangle _\tau} _{\mathcal{V}} & = -  \underbrace{\lim _{\tau \to \infty}  \langle \sum^{N} _{i} m_i \dot{r}_i^2 \rangle _\tau} _{2E_k} \\
   \mathcal{V} = -2\Xi + \mathcal{V}^{ex} & = -2E_k \\
\end{align*}  
$$ 

$\mathcal{V}^{ex}$ 是来自外部对盒子的压力,  $\Xi$  是维里张量

# Compute pressure from virial tensor
$$
\Large
\begin{align*}  
 \mathcal{V}^{ex} & = \langle \sum^{N} _{i} F _{p,i} r_i \rangle _\tau  =  \langle N  \underbrace{ \overline {\sum^{N} _{i} F _{p,i} r_i}  } _{\overline{F _{p,i}} = -p \mathrm{d}S, r=r_i} \rangle _\tau = -p\int r \mathrm{d}S = -3pV \\
 \Xi & = -\sum _{i \lt j} r _{ij} F _{ij} /2\\
  \mathcal{V} & = -2\Xi -3pV  = -2E_k \\
 p & =\frac{2}{3V} (E_k-\Xi)
\end{align*}  
$$ 


# Berendsen method

分子内相互作用对压力的贡献可以忽略, 压力的改变可以转化为对分子间距离的放缩, 同时放缩盒子尺寸. 为此可将运动方程 

$$
\Large
v=\dot{x}
$$ 

修改为

$$
\Large
v + \alpha x= \dot{x} 
$$ 

相应的体积也会变成 

$$
\Large
\dot{V} = 3\alpha V
$$

恒温可压缩系数( $\beta$ )为

$$
\Large
\begin{align*} 
\beta = -\frac{1}{V} \frac{\partial V}{\partial p} \\
\frac{\mathrm{d} p}{\mathrm{d} t} = \dot{p}= -\frac{1}{\beta V} \dot{V} = -\frac{3\alpha}{\beta }
\end{align*}  
$$

带入上式到运动方程:

$$
\Large
\dot{x} = v + \frac{\beta (p_0 - p)}{3 \tau_p}x ; ~~~\dot{p} = (p_0 - p) / \tau_p \\
$$

等式两端同除以 $v$ 可得缩放因子 $\mu$ :

$$
\Large
\mu = 1 + \frac{\beta \Delta t (p_0 - p)}{3 \tau_p} ; \\
$$

上式中 $\Delta t$ 是模拟步长.

