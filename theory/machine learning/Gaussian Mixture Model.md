# Gaussian Mixture Model

### 如何计算高斯混合模型中某个点对应的概率

```python
#scikit-learn
import numpy as np
from sklearn.mixture import GaussianMixture
X = np.array([[1, 2], [1, 4], [1, 0], [10, 2], [10, 4], [10, 0]])
gm = GaussianMixture(n_components=2, random_state=0).fit(X)
gm.means_
gm.predict_proba([[0, 0], [12, 3]])
```
在高斯混合模型拟合数据后，体系会被分成 $c$ 个高斯成分，其中每个高斯成分( $G_i$ )的所占的权重( $w_i$ )各不相同。若体系有 $n$ 个描述符，则某个 $n$ 维空间中点( $r$ )属于 ($G_i$) 的概率 $P_i (r)$ 为：

$$
P_i (r) =  w_i G_i(r, \mu_i, \sigma_i) / \sum^{c}_{i=1} w_i G_i(r, \mu_i, \sigma_i) 
$$

其中 $\mu_i$ , $\sigma_i$ 为 $G_i$ 的均值和协方差矩阵

实际代码运算中均是以对数进行运算的，以三个二维高斯函数为例：

```python
def gaussian_2d(x, y, mean, cov):  
    inv_cov = np.linalg.inv(cov)  
    det_cov = np.linalg.det(cov)  
    dx = x - mean[0]  
    dy = y - mean[1]  
    exponent = dx * inv_cov[0,0] * dx + dy * inv_cov[1,1] * dy + 2 * (dx * inv_cov[0,1] * dy)  
    return 1 / (2 * np.pi * np.sqrt(det_cov)) * np.exp(-0.5 * exponent)  

w=GMM.weights_
mu=GMM.means_
sigma=GMM.covariances_
x=0.37
y=-0.33

#wG is weighted Gaussians
wG= np.log(w)+np.log(
                [gaussian_2d(x,y, mu[0], sigma[0]),
                 gaussian_2d(x,y, mu[1], sigma[1]),
                 gaussian_2d(x,y, mu[2], sigma[2])
                ]
            )

P = np.exp(wG-np.log(np.sum(np.exp(wG))))

```
