#!/bin/bash



gmx_mpi grompp -f ./MDP/prd.mdp -c ../eq.gro -p spce.top -o ./rerun/spce.tpr 

for i in `seq 0 0.1 1`
do 


gmx_mpi mdrun -v -s ./rerun/spce.tpr -deffnm ./rerun/$i -rerun ./md/$i.xtc 

coul=`echo "3 4 0"|gmx_mpi energy -f ./rerun/$i.edr -o ./rerun/$i.xvg |grep Coul|awk '{ total += $3 } END { print total }' `

echo "$coul" >> energy.dat


done
    
