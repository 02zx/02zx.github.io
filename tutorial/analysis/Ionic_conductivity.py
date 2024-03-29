import MDAnalysis as mda
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import sys

fnm=sys.argv[1]
Temp=308.15

#---read trajectories-------------
tpr = '{}.gro'.format(fnm)
traj = '{}.xtc'.format(fnm)
selection_cation = 'name LP'
selection_anion = 'name B'
selection_ion = 'name LP or name B'

u = mda.Universe(tpr, traj)
N=u.atoms.select_atoms(selection_ion).n_atoms
V = u.dimensions[0] * u.dimensions[1] * u.dimensions[2]/1000
total_num_frames = len(u.trajectory)
print('# of frames:', total_num_frames)
print('# of ions:', N)
print('reading gromacs trajectory... ')
print('time (ps), length (nm)')
coords=np.zeros((N,3, total_num_frames),dtype=float)
i=0
for ts in u.trajectory:
    coords[:,:,i] = u.atoms.select_atoms(selection_ion).positions
    i=1+i
    
#--------------------------------------

print('#####################################################')
print('                                                    #')
print('    You must use -pbc nojump to convert trajectory! #')
print('                                                    #')
print('#####################################################')

#----Settings-------------------------------------------------------------------------
unit_t=1 #set the time interval between each frame, ps
T=5000        #ps
dt=1     #ps
ave_time=5000  #ps
print('settings:')
print('computing trajectory in {} ps/frame'.format(unit_t))
print('computing time from 0 to {} ps with a interval of {} ps'.format(T, dt))
print('each point are averaged by {} ps'.format(ave_time))
print('Temperature: {} K'.format(Temp))

ave_fr = int(np.round(ave_time/unit_t))

length = int(np.round(T/dt))

if ave_fr+T/dt >= total_num_frames:
    print("Error! trajectory is not long enough. please reduce T or ave_time")
#--------------------------------------------------------------------------------

t = np.arange(1,T,dt)
dfr = np.round(t/unit_t)

sep=int(N/2)
def conductivity(dfr):     
    #unit nm^2
    upper_triangle = np.zeros((ave_fr,N,N), dtype=float) 
    diagonal_sum = np.zeros((ave_fr,), dtype=float)  
    self_sum = np.zeros((ave_fr,), dtype=float)
    cross_sum = np.zeros((ave_fr,), dtype=float)
    for fr in range(ave_fr):  
            upper_triangle[fr,:,:]=np.triu(np.dot((coords[:,:, fr + int(dfr)] - coords[:,:, fr]),(coords[:,:, fr + int(dfr)] - coords[:,:, fr]).T))
#            upper_triangle = np.triu(np.dot(vecs, vecs.T))
            diagonal_sum[fr] = np.trace(upper_triangle[fr,:,:])
            self_sum[fr] = np.sum(upper_triangle[fr,:sep,:sep])+np.sum(upper_triangle[fr,sep:,sep:])-diagonal_sum[fr]
            cross_sum[fr]=np.sum(upper_triangle[fr,:sep,sep:])
    return np.mean(diagonal_sum)/100, np.mean(self_sum)/100, -1*np.mean(cross_sum)/100



print('computing ionic conductivity...')
with ProcessPoolExecutor() as executor:  
    IC_result = list(tqdm(executor.map(conductivity, dfr), total=length)) 
    
    
ionic_contribution=np.array(IC_result)
ionic_contribution[0,:]


#output

ionic_conductivity=np.sum(ionic_contribution, axis=1)

#fitting

slope_total, _ = np.polyfit(t,ionic_conductivity*(1.86e+6)/6/V/Temp, 1)
slope_NE, _ = np.polyfit(t,ionic_contribution[:,0]*(1.86e+6)/6/V/Temp, 1)

with open('{}.conduct'.format(fnm), 'w') as f:
    print('#Total conductivity, NE conductivity (S/m): {} {}'.format(round(slope_total,3), round(slope_NE,3)), file=f)
    print('#time(ps) total_conductivity(S*ps/m) sum of msd for all N particles (nm^2) self-contribution(nm^2) cross-contribution(nm^2)', file=f)
    print('#Volume: {} nm^3'.format(round(V,3)), file=f)
    print('#Temperature: {} K'.format(round(Temp,3)), file=f)
    print('#contribution(nm^2)--*(1.86e+6)/6/V/Temp-->conductivity(S*ps/m)', file=f)
    print('0 0 0 0 0', file=f)
    for i in range(len(t)):
        print("{:.2f} {:.4g} {:.4g} {:.4g} {:.4g}".format(t[i],ionic_conductivity[i]*(1.86e+6)/6/V/Temp, ionic_contribution[i,0],ionic_contribution[i,1],ionic_contribution[i,2] ), file=f)



