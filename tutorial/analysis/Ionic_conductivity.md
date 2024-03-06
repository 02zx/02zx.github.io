
The ionic conductivity $\sigma$ can be computed from:

$$
\sigma = \frac{1}{6k_BVT} \lim_{t-> \infty} \frac{\rm d}{\rm dt} \bigg{\langle} \sum_{i=1}^{N}  \sum_{j=1}^{N} \[ \vec{r}_i(t) -\vec{r}_i(0) \] \cdot \[ \vec{r}_j(t) - \vec{r}_j(0) \] \bigg{\rangle}
$$

where V is the volume of simulation box, N is the number of ions. This equation takes into consideration all cation−cation,anion−anion, and cation−anion interactions and averagesthose interactions over all frames of the simulation. 
