import MDAnalysis as mda
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import sys

#---read trajectories-------------
tpr = '/home/Share/XZ/conf.pdb'
fnm=sys.argv[1]
traj = '{}.lammpstrj'.format(fnm)
nbins=int(sys.argv[2])

#-----------------------------------------------
u = mda.Universe(tpr, traj,format='LAMMPSDUMP')
selection = 'name O'
N=u.atoms.select_atoms(selection).n_atoms

total_num_frames = len(u.trajectory)
print('# of frames:', total_num_frames)
print('# of particles:', N)
print('reading gromacs trajectory... ')
print('time (ps), length (nm)')

coords=np.zeros((N,3, total_num_frames),dtype=float)
i=0
for ts in u.trajectory:
    coords[:,:,i] = u.atoms.select_atoms(selection).positions
    i=1+i
    
    
#--------------------------------------

print('#####################################################')
print('                                                    #')
print('    You must use xu,yu,zu for dump settings!        #')
print('                                                    #')
print('#####################################################')
#----input-------------------------------------------------------------------------
unit_t=0.1 #set the time interval between each frame, ps
T=5        #ps
dt=0.1     #ps
ave_time=5  #ps
print('settings:')
print('computing trajectory in {} ps/frame'.format(unit_t))
print('computing time from 0 to {} ps with a interval of {} ps'.format(T, dt))
print('each point are averaged by {} ps'.format(ave_time))


ave_fr = int(np.round(ave_time/unit_t))

length = int(np.round(T/dt))

if ave_fr+T/dt >= total_num_frames:
    print("Error! trajectory is not long enough. please reduce T or ave_time")
#--------------------------------------------------------------------------------

t = np.arange(1*dt,T,dt)
dfr = np.round(t/unit_t)
bin_z=np.max(coords[:,2,:])/nbins
def mean_square_displacement(dfr):  
    result = np.zeros((ave_fr,nbins), dtype=float)  
    for fr in range(ave_fr):  #loop for time
        for z in range(nbins):#loop for z position, take sum for each bin
            result[fr,z]=np.sum(np.square(coords[:,:, fr + int(dfr)] - coords[:,:, fr])[(coords[:,2,fr]>=0+z*bin_z) & (coords[:,2,fr]<0+(z+1)*bin_z)])/N           
    return np.mean(result,axis=0)  #take average for time


print('computing mean square displacement...')
with ProcessPoolExecutor() as executor:  
    msd = list(tqdm(executor.map(mean_square_displacement, dfr), total=length))  

data=np.array(msd)
output = np.zeros((nbins,2), dtype=float)
for num in range(nbins):
    output[num,1], _ = np.polyfit(t,data[:,num], 1) #slope ->second col of output
    output[num,0]=0+num*bin_z+bin_z/2 # z coordinate -> first col of output

with open('{}.dz'.format(fnm), 'w') as f:
    print('#z (nm), Diffusivity (nm^2/ps)', file=f)
    for i in range(nbins):
        print("{:.3f} {:.4g}".format(output[i,0]/10, output[i,1]/100), file=f)

with open('msd.log','w') as f:
    print('#time(ps) msd(bin1->nbins) (nm^2)', file=f)
    for i in range(len(t)):
        print("{:.3f}".format(t[i]), end=' ',file=f)
        for j in range(nbins):
            print("{:.4g}".format(data[i,j]), end=' ',file=f)
        print(" ",file=f)
