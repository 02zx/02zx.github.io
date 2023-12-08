#!/bin/bash
mkdir em eq md rerun
for i in `seq 0 0.1 1` #lambda = 0, 0.1, 0.2, ... 1
do 
#compute charge
hw=`echo "$i"|awk '{print 0.4238*sqrt($1)}'` 
ow=`echo "$i"|awk '{print -2*0.4238*sqrt($1)}'`

#modify top file
sed -i s"@OW_spc.*@OW_spc       8      15.9994  $ow  A   3.16557e-01  6.50629e-01@" topo.top
sed -i s"@HW_spc.*@HW_spc       1       1.0080  $hw  A   0.00000e+00  0.00000e+00@" topo.top

#em                                                                    generate new top to em/
gmx_mpi grompp -f MDP/em.mdp -c ../eq.gro -p topo.top -o em/$i -maxwarn 3 -pp em/$i.top
gmx_mpi mdrun -v -deffnm em/$i 

sed -i s"@i=.*@i=$i@" run.sh
#eq/md
sbatch run.sh


done

    
