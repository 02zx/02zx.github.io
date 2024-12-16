## 模拟超声状态

### 源代码修改记录：
本代码基于gromacs2023.5版本进行修改，将Berendsen的半各向异性控压(semiisotropic)改为了各向同性控压(isotropic)+超声波(压力的正弦变化)。

压力变化遵循以下公式：

$$ 
p(t) = p_0 + A*\mathrm{sin}(\frac{2\pi t}{T}) 
$$

其中 $p(t)$ 是t时刻的压力， $p_0$ 是无超声时的参考压力， $A$ 是压力的振幅， $T$ 是压力振动周期 
### 使用方法：

从官网下载gromacs-2023.5并解压。将[coupling.cpp](https://github.com/02zx/02zx.github.io/blob/main/tutorial/molecular%20dynamics%20simulation/coupling.cpp)复制到gromacs-2023.5/src/gromacs/mdlib/coupling.cpp之后进行编译即可。


mdp设置：

mdp中的参数与以上公式相对应，其中 $T$ 是以步数的形式确定的，例如 $T=500000$ 时， 在2fs步长下，振动周期即为 500000*2fs = 1ns.
```
Pcoupl              =  Berendsen ; for ultrasonic simulation
Pcoupltype          =  semiisotropic ;
tau_p               =  1
compressibility     =  4.5e-5 500000 ; compressibility,  T(steps)
ref_p= 0 100; p_0(bar), A(bar)
```


