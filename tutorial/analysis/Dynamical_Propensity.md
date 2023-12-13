# Dynamical Propensity of Liquid Water

https://doi.org/10.1073/pnas.1817135116

1. 进行NVT模拟, 计算O-O的径向分布函数 $g_{OO}$ (vmd或gmx计算获得)

2. 计算isotropic structure factor $S(q)$ , 并确定第一个峰出现的位置 $q_0$ 


$$
S(q) = 1+ \rho \int_{0}^{\infty}  e^{-i\bf{k}r} g_{OO}(r) dr
$$

$$
S(q) = 1+ \frac{4\pi \rho}{q} \int_{0}^{\infty}  r \sin(rq) [g_{OO}(r)-1] dr
$$

3. 计算 $\Phi (\bf{q},t)$

$$
 \Phi (\bf{q},t) = \frac{1}{N} \sum_{j=1}^{N} e^{i\bf{q} [r_j(t) - r_j(0)]}
$$

4. 计算 self-intermediate scattering function $F(\bf{q},t)$ 和 dynamical susceptibility $\chi_4 (\bf{q},t)$ 计算200个各向同性方向上的平均值

$$
F(\bf{q},t) = \langle  \Phi (\bf{q},t) \rangle _{|q| = q_0}
$$

$$
\chi_4 (\bf{q},t) = N \Big{[} \langle | \Phi (\bf{q},t) |^2  \rangle - \langle\Phi (\bf{q},t)\rangle^2   \Big{]} 
$$

5. The time of maximum heterogeneity $t_0$ 取自 $\chi_4 (q_0,t)$ 最大值出现的位置

6. 计算dynamical propensity

$$
DP_i = \Big{\langle}   \frac{|\bf{r}_i(t_0) - \bf{r}_i(0)|^2}{MSD}   \Big{\rangle} _{ISO}
$$



# Appendix
## For a periodic function f(t) with periodicity of T.

$$
f(t) = A_{0}+\sum_{n=1}^{\infty} [A_n cos(2\pi n t /T) + B_n sin(2\pi n t /T)]
$$

### Orthogonal functions

 $\{sin(0x),cos(0x),sin(x),cos(x),sin(2x),cos(2x),...,sin(nx),cos(nx)\}$

$$
\int_{-\pi}^{\pi} sin(nx) cos(mx) = 0
$$

For n $\ne$ m

$$
\int_{-\pi}^{\pi} cos(nx) cos(mx) = 0
$$

$$
\int_{-\pi}^{\pi} sin(nx) sin(mx) = 0
$$



###  Fourier coefficient 

$x=2\pi t /T$ 

$$
f(x) = A_{0}+ \sum_{n=1}^{\infty} [A_n cos(nx) + B_n sin( nx)]
$$

### #1 For n = 0


$$
\int_{-\pi}^{\pi} f(x) cos(0x) dx = \int_{-\pi}^{\pi} [ A_{0}+\sum_{n=1}^{\infty} (A_n cos( n x ) + B_n sin( n x )) ] cos(0x) dx = \int_{-\pi}^{\pi}( A_{0} cos(0x) + 0 + 0 )dx = \int_{-\pi}^{\pi} A_{0} dx = A_{0}2\pi
$$

$$
\Rightarrow A_0 = (1/2\pi) \int_{-\pi}^{\pi} f(x)dx = (1/T) \int_{-T/2}^{T/2} f(t)dt
$$

### #2 For n $\ge$ 1


$$
\int_{-\pi}^{\pi} f(x) cos(nx) dx = \int_{-\pi}^{\pi} [ A_{0}+\sum_{n=1}^{\infty} (A_n cos( n x ) + B_n sin( n x )) ] cos(nx) dx = \int_{-\pi}^{\pi}( 0+ A_n cos( n x ) cos(nx) + 0 )dx = \int_{-\pi}^{\pi} A_n cos^2(nx) dx = A_n \pi
$$

$$
\Rightarrow A_n = (1/\pi) \int_{-\pi}^{\pi} f(x)cos(nx) dx = (2/T) \int_{-T/2}^{T/2} f(t) cos(2\pi n t /T) dt 
$$

### #3 For n $\ge$ 1
$$
\int_{-\pi}^{\pi} f(x) sin(nx) dx = \int_{-\pi}^{\pi} [ A_{0}+\sum_{n=1}^{\infty} (A_n cos( n x ) + B_n sin( n x )) ] sin(nx) dx = \int_{-\pi}^{\pi}( 0+ 0+ B_n sin(nx) sin(nm) )dx = \int_{-\pi}^{\pi} B_n sin^2(nx) dx = B_n \pi
$$

$$
\Rightarrow B_n = (1/\pi) \int_{-\pi}^{\pi} f(x) sin(nx) dx= (2/T) \int_{-T/2}^{T/2} f(t) sin(2\pi n t /T) dt 
$$


### Exponentials

$cos(x) = (e^{ix} + e^{-ix})/2$

$sin(x) = -(e^{ix} - e^{-ix})/2$

$$
f(t) = A_{0}+\sum_{n=1}^{\infty} [A_n cos(2\pi n t /T) + B_n sin(2\pi n t /T)] \Rightarrow f(t) = \sum_{n=-\infty}^{\infty} [(1/T) \int^{T}_{0} f(t) e^{-i2\pi nt/T} dt] e^{i2\pi nt/T}
$$

### Fourier transform F(n) 
$$
F(n) = \int^{T}_{0} f(t) e^{-i2\pi nt/T} dt
$$

### For a discrete function
$$
\{a_m \} = \{a_0, a_1, a_2,..., a_{N-1}\}
$$

$$
t = m \times dt;~ f(t) = \{a_m \};~ T = N \times dt;~~ n = k/T ; \int^{T}_{0} dt \Rightarrow \sum^{N-1}_{m=0}
$$

$$
F(k) = \sum^{N-1}_{m=0} a_m e^{-2\pi i mk/N}~,k=0,1,2,...,N-1
$$
