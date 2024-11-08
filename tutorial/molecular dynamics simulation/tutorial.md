
## NVE Simulation


```python
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import distance

#molecular mechanics
def LJ_interactions(sigma, epsilon, coord):
  N=len(coord)
  rij=coord[:,np.newaxis,:]-coord[np.newaxis,:,:]
  r=distance.squareform(distance.pdist(coord))
  attraction = (sigma/r)**6
  replusion = attraction**2
  potential = 4*epsilon*(replusion - attraction)
  #F=-dV/dr
  force = np.multiply((4*epsilon*(12*replusion/r - 6*attraction/r))[:,:,np.newaxis], np.divide(rij,r[:, :, np.newaxis]))
  return potential, force

#molecular dynamics: Velocity Verlet
def md_vv(sigma,epsilon,mass,coord,velocity,dt):
  N=len(coord)
  force=LJ_interactions(sigma, epsilon, coord)[1]
  #third law
  accelerate=np.zeros((N,3),dtype=float)
  for i in range(N): 
    accelerate[i,:]= np.sum(np.delete(force[i,:,:],i,axis=0),axis=0)/mass
  new_coord=coord + velocity*dt + 0.5*accelerate*dt**2
  new_potential, new_force=LJ_interactions(sigma, epsilon, new_coord)
  new_acc = np.zeros((N,3),dtype=float)
  for i in range(N): 
    new_acc[i,:]= np.sum(np.delete(new_force[i,:,:],i,axis=0),axis=0)/mass
  new_vel = velocity + 0.5*(accelerate+new_acc)*dt
  Ek=np.sum(mass*new_vel**2)/2
  return new_coord, new_vel,Ek,np.sum(np.triu(new_potential,k=1))
```


```python
coord=np.array([[0,1,0],[0,2,0],[0,-0.5,1.8],[0,1,-1.2]])
velocity=np.array([[0,0.001,0],[0,0.002,0],[0,0,0.001],[0,0,0]])
sigma=1.0
epsilon=10000
dt=0.0005
mass=20
kinetic=[]
potential=[]
N=len(coord)

nsteps=8000
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
