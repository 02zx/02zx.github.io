import MDAnalysis as mda
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

#---read trajectories-------------
tpr = 'ow.gro'
traj = 'ow.xtc'
u = mda.Universe(tpr, traj)
N=u.atoms.n_atoms  
V = u.dimensions[0] * u.dimensions[1] * u.dimensions[2]
total_num_frames = len(u.trajectory)
print('# of frames:', total_num_frames)
print('# of particles:', N)
print('reading gromacs trajectory... ')
print('time (ps), length (nm)')

coords=np.zeros((N,3, total_num_frames),dtype=float)
i=0
for ts in u.trajectory:
    coords[:,:,i] = u.atoms.positions
    i=1+i
#--------------------------------------
print('#########################################################')
print('#                                                       #')
print('#  You should use -pbc nojump to convert trajectory!!!  #')
print('#                                                       #')
print('#########################################################')
#----input-------------------------------------------------------------------------
unit_t=0.01 #set the time interval between each frame, ps
T=20        #ps
dt=0.01     #ps
ave_time=20  #ps
print('settings:')
print('computing trajectory in {} ps/frame'.format(unit_t))
print('computing time from 0 to {} ps with a interval of {} ps'.format(T, dt))
print('each point are averaged by {} ps'.format(ave_time))


ave_fr = int(np.round(ave_time/unit_t))

length = int(np.round(T/dt))

if ave_fr+T/dt >= total_num_frames:
    print("Error! trajectory is not long enough. please reduce T or ave_time")
#--------------------------------------------------------------------------------

t = np.arange(1,T,dt)
dfr = np.round(t/unit_t)

def mean_square_displacement(dfr):  
    result = np.zeros((ave_fr,), dtype=float)  
    for fr in range(ave_fr):  
            result[fr] = np.sum(np.square(coords[:,:, fr + int(dfr)] - coords[:,:, fr]))/N
    return np.mean(result)  

def mean_fouth_order_displacement(dfr):  
    result = np.zeros((ave_fr,), dtype=float)  
    for fr in range(ave_fr):  
            result[fr] = np.sum(np.power((coords[:,:, fr + int(dfr)] - coords[:,:, fr]),4))/N
    return np.mean(result)  

print('computing mean square displacement...')
with ProcessPoolExecutor() as executor:  
    msd = list(tqdm(executor.map(mean_square_displacement, dfr), total=length))  

    
print('computing mean fouth order displacement...')

with ProcessPoolExecutor() as executor:  
    mfd = list(tqdm(executor.map(mean_fouth_order_displacement, dfr), total=length))  


def alpha_2(msd,mfd):
    return 3*np.array(mfd)/5/(np.array(msd)**2)-1
 

print('computing alpha_2...')

a_2 = alpha_2(msd, mfd)
with open('alpha_2.txt', 'w') as f:
    print('0 0 0 0', file=f)
    for i in range(len(t)):
        print("{:.2f} {:.4g} {:.4g} {:.4g}".format(t[i], a_2[i], msd[i]/100, mfd[i]/10000), file=f)
