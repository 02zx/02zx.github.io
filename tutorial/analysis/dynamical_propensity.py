import numpy as np
import scipy.spatial as spatial
import networkx as nx
import matplotlib.pyplot as plt
import ase
from ase.io import read
from MDAnalysis.analysis.distances import distance_array
import MDAnalysis as mda
from MDAnalysis.coordinates import XTC
import sys

MSD=0.117553*100
ref=mda.Universe('initial.pdb')
ref_atm=ref.atoms
ref_coords=ref_atm.positions
box_size=ref.dimensions
dr=np.zeros((1,13824))
for i in range(1,6):
    fnm='{}/md.pdb'.format(i)
#fnm='../0/OW'


    u=mda.Universe('{}'.format(fnm))

    atoms=u.atoms
    coords=atoms.positions
    delta=coords-ref_coords
    delta[:,0]=delta[:,0]-box_size[0]*np.round(delta[:,0]/box_size[0])
    delta[:,1]=delta[:,1]-box_size[0]*np.round(delta[:,1]/box_size[0])
    delta[:,2]=delta[:,2]-box_size[0]*np.round(delta[:,2]/box_size[0])
    
        
    dr=dr+np.square(np.linalg.norm(delta,axis=1))
    
DH=dr/4/MSD

print('%5:',np.sort(DH[0])[692], '%95', np.sort(DH[0])[13133])

with open('DH','w') as f:
    for i in DH[0]:
        print(round(i,3),file=f)
