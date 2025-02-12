
# Molecular Dynamics Simulation

## NVE Ensemble
### Equation of Motions
#### Velocity Verlet
$r(t+\delta t) = r(t) + v(t)\delta t + 0.5a(t) \delta t^2$

$v(t+\delta t) = v(t) +  0.5(a(t) +a(t+\delta t))\delta t$

```python
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import distance

#molecular mechanics: Lennard-Jones potential
def LJ_interactions(sigma, epsilon, coord):
  N=len(coord)
  rij=coord[:,np.newaxis,:]-coord[np.newaxis,:,:] #NxNx3 matrix of coordinates
  r=distance.squareform(distance.pdist(coord)) #NxN distance matrix
  attraction = (sigma/r)**6 #NxN
  replusion = attraction**2 #NxN
  potential = 4*epsilon*(replusion - attraction) #NxN
  #F=-dV/dr
  force = np.multiply((4*epsilon*(12*replusion/r - 6*attraction/r))[:,:,np.newaxis], np.divide(rij,r[:, :, np.newaxis])) #NxNx3 
  return potential, force

#molecular dynamics: Velocity Verlet
def md_vv(sigma,epsilon,mass,coord,velocity,dt):
  N=len(coord)
  force=LJ_interactions(sigma, epsilon, coord)[1]

  #Second law:F=ma----
  accelerate=np.zeros((N,3),dtype=float)
  for i in range(N): 
    accelerate[i,:]= np.sum(np.delete(force[i,:,:],i,axis=0),axis=0)/mass #Nx3 acceleration for each particle 
  
  #update coordinates & velocity----
  new_coord=coord + velocity*dt + 0.5*accelerate*dt**2 #Nx3
  new_potential, new_force=LJ_interactions(sigma, epsilon, new_coord) 
  new_acc = np.zeros((N,3),dtype=float)
  for i in range(N): 
    new_acc[i,:]= np.sum(np.delete(new_force[i,:,:],i,axis=0),axis=0)/mass
  new_vel = velocity + 0.5*(accelerate+new_acc)*dt

  #Compute Kinetic energy-----------
  Ek=np.sum(mass*new_vel**2)/2 
  return new_coord, new_vel,Ek,np.sum(np.triu(new_potential,k=1))
```


```python
#initial configurations
coord=np.array([[0,1,0],[0,2,0],[0,-0.5,1.8],[0,1,-1.2]])
velocity=np.array([[0,0.001,0],[0,0.002,0],[0,0,0.001],[0,0,0]])

#forcefield
sigma=1.0
epsilon=10000
mass=20

#integration 
dt=0.0005


kinetic=[]
potential=[]
N=len(coord)


traj=np.zeros((nsteps,N,3),dtype=float)
for i in range(nsteps):
  coord, velocity, k,pe=md_vv(sigma,epsilon,mass,coord,velocity,dt)
  kinetic.append(k)
  potential.append(pe)
  traj[i,:,:]=coord
  print(i,k,pe)

plt.plot(range(nsteps),potential)
plt.plot(range(nsteps),kinetic)

plt.plot(traj[:-1,1,1],traj[:-1,1,2],'-')
plt.plot(traj[:-1,2,1],traj[:-1,2,2],'-')
plt.plot(traj[:-1,0,1],traj[:-1,0,2],'-')
plt.plot(traj[:-1,3,1],traj[:-1,3,2],'-')
```

#### Periodic Boundary Condition

$\delta r= r_2 - r_1 = (\delta x, \delta y, \delta z)$

$\delta r_{pbc} = \delta r - cell* round(\delta r/cell)$

$cell$ 是一个 3 $\times$ 3 矩阵, 其非对角项均为0, 对角项分别为盒子的xyz边长。 目前该pbc程序只支持矩形盒子。

```python
def pbc(coord,cell):
  rij=coord[:,np.newaxis,:]-coord[np.newaxis,:,:]
  rij_pbc=rij-np.diag(cell)[np.newaxis,np.newaxis,:]*np.round(np.divide(rij,np.diag(cell)[np.newaxis,np.newaxis,:]))
  r_pbc=np.linalg.norm(rij, ord=2,axis=2)
  return rij_pbc, r_pbc

def LJ_interactions_pbc(sigma, epsilon, coord, cell):
  N=len(coord)
  #--use PBC--------------
  rij, r=pbc(coord,cell) 
  attraction = (sigma/r)**6
  replusion = attraction**2
  potential = 4*epsilon*(replusion - attraction)
  #F=-dV/dr
  force = np.multiply((4*epsilon*(12*replusion/r - 6*attraction/r))[:,:,np.newaxis], np.divide(rij,r[:, :, np.newaxis]))
  return potential, force

def md_vv_pbc(sigma,epsilon,mass,cell,coord,velocity,dt):
  N=len(coord)
  #--use PBC--------------
  force=LJ_interactions_pbc(sigma, epsilon, coord)[1]
  #third law
  accelerate=np.zeros((N,3),dtype=float)
  for i in range(N): 
    accelerate[i,:]= np.sum(np.delete(force[i,:,:],i,axis=0),axis=0)/mass
  new_coord=coord + velocity*dt + 0.5*accelerate*dt**2
  #--use PBC--------------
  new_potential, new_force=LJ_interactions_pbc(sigma, epsilon, new_coord)
  new_acc = np.zeros((N,3),dtype=float)
  for i in range(N): 
    new_acc[i,:]= np.sum(np.delete(new_force[i,:,:],i,axis=0),axis=0)/mass
  new_vel = velocity + 0.5*(accelerate+new_acc)*dt
  Ek=np.sum(mass*new_vel**2)/2
  return new_coord, new_vel,Ek,np.sum(np.triu(new_potential,k=1))
```

```python
#initial configurations
xx=10;yy=10;zz=10
xy=0;xz=0;yx=0;yz=0;zx=0;zy=0
boxsize=np.array([
[xx,xy,xz],
[yx,yy,yz],
[zx,zy,zz],
])

```

#### Coulombic Interactions

### Calculation of Properties
#### Temperature
#### Pressure
#### Virial
#### Diffusion



## NVT Ensemble
### Berendsen Thermostat
### Velocity Rescale Thermostat
### Anderson Thermostat
### Nose-Hoover Thermostat

## NPT Ensemble
### Berendsen Barostat
### Parrinello-Rahman Barostat
### C-rescale Barostat

