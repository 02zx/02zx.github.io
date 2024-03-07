
The ionic conductivity $\sigma$ can be computed from:

$$
\sigma = \frac{1}{6k_BVT} \lim_{t-> \infty} \frac{\rm d}{\rm dt} \bigg{\langle} \sum_{i=1}^{N}  \sum_{j=i}^{N}  q_iq_j\[ \vec{r}_i(t) -\vec{r}_i(0) \] \cdot \[ \vec{r}_j(t) - \vec{r}_j(0) \] \bigg{\rangle}
$$

where V is the volume of simulation box, N is the number of ions. This equation takes into consideration all cation−cation,anion−anion, and cation−anion interactions and averagesthose interactions over all frames of the simulation.
Above equation can be simplified as:

$$
\sigma = \sigma_{NE} + \sigma_{self} + \sigma_{cross}
$$

$$
\sigma_{NE} = \frac{N_{pair}}{k_BVT} \bigg{(} \sum_i^{k} q_i^2 \nu_i D_i \bigg{)} 
$$

$$
\sigma_{self} = \frac{1}{6k_BVT} \sum_{i=1}^k \lim_{t-> \infty} \frac{\rm d}{\rm dt} \bigg{\langle}   \sum_{a=1}^{N_i-1}  \sum_{b=a+1}^{N_i}  q_i^2\[ \vec{r}_a(t) -\vec{r}_a(0) \] \cdot \[ \vec{r}_b(t) - \vec{r}_b(0) \]  \bigg{\rangle} 
$$

$$
\sigma_{cross} = \frac{1}{6k_BVT} \sum_{i=1}^{k-1} \sum_{j=i+1}^k \lim_{t-> \infty} \frac{\rm d}{\rm dt} \bigg{\langle}   \sum_{a=1}^{N_i}  \sum_{b=1}^{N_j}  q_iq_j\[ \vec{r}_a(t) -\vec{r}_a(0) \] \cdot \[ \vec{r}_b(t) - \vec{r}_b(0) \]  \bigg{\rangle} 
$$

where $k$ is the number of ion species, $\nu_i$, $q_i$, $D_i$ are the stoichiometric number, charge, and self-diffusion coefficient for $i$-th specie. 

The first term $\sigma_{NE}$ is the Nernst−Einstein equation (law of the independent migration of ions). 

The second term $\sigma_{self}$ considers all the ion-ion interaction between the same specie, for example, $Na^+ - Na^+$,  $K^+ - K^+$, and $Cl^- - Cl^-$ are considered, but $Na^+ - Na^+$, $K^+ - Cl^-$ and $Na^+ - Cl^-$ are not. 

The third term $\sigma_{cross}$ considers all the ion-ion interaction between the different species, for example, $Na^+ - K^+$, $K^+ - Cl^-$ and $Na^+ - Cl^-$.


 
