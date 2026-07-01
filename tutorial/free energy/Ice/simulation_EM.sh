#2026/6/26

#--Settings-------------------------
Pressure=("0")    #bar
Temperature=("200")  #K


#-------------------
Equilibration=no
#-------------------
Einstein=no
#-------------------
Decoupl=no
    num_points=16
#-------------------
Analyze_Einstein=yes
Analyze_Decoupl=yes




#----------------------------------------------
let length=${#Pressure[@]}-1
lambda0=`echo "1662800"|awk '{print $1/200/2}'` #5000kT/nm^2 T=1K
#equilibration
if [ "${Equilibration}" = "yes" ]; then

mkdir mdp
mkdir src

cat <<EOF >mdp/em.mdp
integrator = steep
nsteps = 10000
emtol  = 10.0
emstep = 0.01
;
nstxout   = 10000
nstlog    = 50
nstenergy = 50
;
pbc = xyz
cutoff-scheme            = Verlet
coulombtype              = PME
rcoulomb                 = 0.85
vdwtype                  = Cut-off
rvdw                     = 0.85
rlist			         = 0.85
;
constraints              = hbonds
;freezegrps          = MOL
;freezedim           = Y Y Y 
;


EOF


cat <<EOF >mdp/eq.mdp

;VARIOUS PREPROCESSING OPTIONS
title                    = equilibration

; RUN CONTROL PARAMETERS
integrator               = md
; Start time and timestep in ps
tinit                    = 0
dt                       = 0.002
nsteps                   = 2500000
; For exact run continuation or redoing part of a run
init_step                = 0
; mode for center of mass motion removal
comm-mode                = linear
; number of steps for center of mass motion removal
nstcomm                  = 100
; group(s) for center of mass motion removal
comm-grps                = system

; OUTPUT CONTROL OPTIONS
; Output frequency for coords (x), velocities (v) and forces (f)
nstxout                  = 5000
nstvout                  = 5000
nstfout                  = 5000
; Output frequency for energies to log file and energy file
nstlog                   = 10000 
nstenergy                = 10000
; Output frequency and precision for xtc file
nstxtcout                = 10000
xtc-precision            = 10000

; NEIGHBORSEARCHING PARAMETERS
; nblist update frequency
nstlist                  = 1
; ns algorithm (simple or grid)
ns_type                  = grid
; Periodic boundary conditions: xyz (default), no (vacuum)
; or full (infinite systems only)
pbc                      = xyz
; nblist cut-off        
rlist                    = 0.85 

; OPTIONS FOR ELECTROSTATICS AND VDW
; Method for doing electrostatics
coulombtype              = pme
rcoulomb                 = 0.85
pme_order                = 4
fourierspacing           = 0.1

; Method for doing Van der Waals
vdw-type                 = cut-off  
rvdw                     = 0.85
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                 = EnerPres 

; OPTIONS FOR WEAK COUPLING ALGORITHMS
; Temperature coupling  
Tcoupl                   = v-rescale
tau_t                    = 1.0
ref_t=200
tc-grps                  = system
Pcoupl              =  Berendsen;c-rescale
Pcoupltype          =  anisotropic
tau_p               =  4       
compressibility     =  4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5 
ref_p=0 0 0 0 0 0 ;off-diagonal don't affact much
; GENERATE VELOCITIES FOR STARTUP RUN
 gen_vel                  = yes
 gen_temp=200
 gen_seed                 = 50934891

; OPTIONS FOR BONDS
constraints              = hbonds ;all-angles
; Type of constraint algorithm
constraint-algorithm     = lincs
lincs-iter               =  4
lincs-order              =  6

EOF


cat <<EOF >topo.top
Include forcefield parameters
#include "amber99sb.ff/forcefield.itp"

;
; Note the strange order of atoms to make it faster in gromacs.
[ atomtypes ]
;name  bond_type    mass    charge   ptype          sigma      epsilon
; tip4pice
  H_ice     1       1.008   0.0000  A   0.00000e+00  0.00000e+00
  O_ice     8       16.00   0.0000  A   3.16680e-01  8.82164e-01
  M     0      0.0000   0.0000  D   0.00000e+00  0.00000e+00
;
[ moleculetype ]
; molname	nrexcl
SOL		2

[ atoms ]
; id	at type	res nr 	residu name	at name	cg nr	charge
1       O_ice        1       SOL      OW     1       0.00
2       H_ice        1       SOL      HW1     1       0.5897
3       H_ice        1       SOL      HW2     1       0.5897
4       M            1       SOL      MW      1       -1.1794

[ settles ]
;i j funct doh  dhh
1     1   0.09572 0.15139

[ exclusions ]
1	2	3	4
2	1	3	4
3	1	2	4
4	1	2	3

; The position of the virtual site is computed as follows:
;
;		O
;  	      
;	    	D
;	  
;	H		H
;
; const = distance (OD) / [ cos (angle(DOH)) 	* distance (OH) ]
;	  0.01577 nm	/ [ cos (52.26 deg)	* 0.09572 nm	]

; Vsite pos x4 = x1 + a*(x2-x1) + b*(x3-x1)

[ virtual_sites3 ]
; Vsite from			funct	a		b
4	1	2	3	1	0.13458335      0.13458335
[ system ]
; Name
ICE in water

[ molecules ]
; Compound        #mols
SOL 1088

EOF



    for i in `seq 0 1 $length`
    do
    nmol=`grep OW eq.gro|wc -l`
    mkdir ${Temperature[${i}]}_${Pressure[${i}]}
    mkdir ${Temperature[${i}]}_${Pressure[${i}]}/eq
    sed '$d' ./topo.top >${Temperature[${i}]}_${Pressure[${i}]}/eq/topo.top
    echo "SOL $nmol" >> ${Temperature[${i}]}_${Pressure[${i}]}/eq/topo.top
    sed -i s"@ref_t=.*@ref_t=${Temperature[${i}]}@" mdp/eq.mdp
    sed -i s"@gen_temp=.*@gen_temp=${Temperature[${i}]}@" mdp/eq.mdp
    sed -i s"@ref_p=.*@ref_p=${Pressure[${i}]} ${Pressure[${i}]} ${Pressure[${i}]} 0 0 0@" mdp/eq.mdp
    gmx grompp -f mdp/em.mdp -c conf.gro -p ${Temperature[${i}]}_${Pressure[${i}]}/eq/topo.top -o ${Temperature[${i}]}_${Pressure[${i}]}/eq/em
    gmx mdrun -v -deffnm ${Temperature[${i}]}_${Pressure[${i}]}/eq/em
    gmx grompp -f mdp/eq.mdp -c em.gro -p ${Temperature[${i}]}_${Pressure[${i}]}/eq/topo.top -o ${Temperature[${i}]}_${Pressure[${i}]}/eq/eq




cat <<EOF >eq.sh 
#!/bin/bash
#SBATCH -n 16
#SBATCH -J eq_${Temperature[${i}]}_${Pressure[${i}]}
#SBATCH -N 1
#SBATCH -p NVIDIAGeForceRTX3090
#SBATCH --gres=gpu:1
#SBATCH -e error.err1
#SBATCH -o output.out1

module purge
module load compiler/gcc/7.3.1
module load compiler/intel/2021.3.0
module load mpi/intelmpi/2021.3.0
module load mathlib/fftw/intelmpi/3.3.9_single
module load apps/gromacs/intelmpi/2018.8-30
gmx mdrun -ntmpi 1 -ntomp 16 -pin on -deffnm ./${Temperature[${i}]}_${Pressure[${i}]}/eq/eq -v -gpu_id 0 -pme gpu -nb gpu
EOF

sbatch eq.sh
    done
fi

#Einstein crystal
if [ "${Einstein}" = "yes" ]; then

cat <<EOF >einstein.mdp
;VARIOUS PREPROCESSING OPTIONS
title                    = EM iceIh

; RUN CONTROL PARAMETERS
integrator               = md
; Start time and timestep in ps
tinit                    = 0
dt                       = 0.001
nsteps                   = 5000000
; For exact run continuation or redoing part of a run
init_step                = 0
; mode for center of mass motion removal
comm-mode                = none
; number of steps for center of mass motion removal
nstcomm                  = 1
; group(s) for center of mass motion removal
comm-grps                =

; OUTPUT CONTROL OPTIONS
; Output frequency for coords (x), velocities (v) and forces (f)
nstxout                  = 5000
nstvout                  = 5000
nstfout                  = 5000
; Output frequency for energies to log file and energy file
nstlog                   = 10000 
nstenergy                = 10000
; Output frequency and precision for xtc file
nstxtcout                = 10000
xtc-precision            = 10000

; NEIGHBORSEARCHING PARAMETERS
; nblist update frequency
nstlist                  = 1
; ns algorithm (simple or grid)
ns_type                  = grid
; Periodic boundary conditions: xyz (default), no (vacuum)
; or full (infinite systems only)
pbc                      = xyz
; nblist cut-off        
rlist                    = 0.85 

; OPTIONS FOR ELECTROSTATICS AND VDW
; Method for doing electrostatics
coulombtype              = pme
rcoulomb                 = 0.85
pme_order                = 4
fourierspacing           = 0.1

; Method for doing Van der Waals
vdw-type                 = cut-off  
rvdw                     = 0.85
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                 = EnerPres 

; OPTIONS FOR WEAK COUPLING ALGORITHMS
; Temperature coupling  
Tcoupl                   = v-rescale
tau_t                    = 0.2
ref_t                    = 250.
tc-grps                  = system

Pcoupl                   =  no 

; GENERATE VELOCITIES FOR STARTUP RUN
 gen_vel                  = yes
 gen_temp                 = 250.00
 gen_seed                 = -1

; OPTIONS FOR BONDS
constraints              = all-angles
; Type of constraint algorithm
constraint-algorithm     = lincs
lincs-iter               =  6
lincs-order              =  8

freezegrps = SWF1
freezedim  = y y y 
EOF

cat <<EOF >src/ideal.template
[ defaults ]
; nbfunc    comb-rule   gen-pairs   fudgeLJ fudgeQQ
  1     2       no      1.0 1.0

[atomtypes]
;name     mass      charge   ptype    sigma        epsilon
MW     0             0.000       D   0.00000       0.00000
OW     15.99940      0.000       A   0.00000       0.00000
HW     1.00800       0.000       A   0.00000       0.00000
SW     15.99940      0.000       A   0.00000       0.00000

[ moleculetype ]
; molname       nrexcl
S                1

[atoms]
; nr type resnr residu atom cgnr charge
1     SW   1     S  SWF1  1     0         15.994
2     HW   1     S  HWF2  1     0.0000    1.008
3     HW   1     S  HWF3  1     0.0000    1.008
4     MW   1     S  MWF4  1     0.0000    0.0

[ position_restraints ]
; ai  funct  fcx    fcy    fcz
   2    1    spring_const. spring_const. spring_const. ; restrains to a point
   3    1    spring_const. spring_const. spring_const. ; restrains to a point

[constraints]
;i j funct doh  dhh
1       2       1       0.09572
1       3       1       0.09572
2       3       1       0.15139

[exclusions]
1       2       3       4
2       1       2       3
3       1       2       4
4       1       2       3

[dummies3]
; Dummy from            funct   a       b
4	1	2	3	1	0.13458335      0.13458335

[moleculetype]
; name nrexcl
water  1

[atoms]
; nr type resnr residu atom cgnr charge
1     OW   1     water  OW1  1     0         15.994 
2     HW   1     water  HW2  1     0.0000    1.008
3     HW   1     water  HW3  1     0.0000    1.008
4     MW   1     water  MW4  1     0.0000    0.0

[ position_restraints ]
; ai  funct  fcx    fcy    fcz
   1    1    spring_const. spring_const. spring_const. ; restrains to a point
   2    1    spring_const. spring_const. spring_const. ; restrains to a point
   3    1    spring_const. spring_const. spring_const. ; restrains to a point

[constraints]
;i j funct doh  dhh
1       2       1       0.09572
1       3       1       0.09572
2       3       1       0.15139

[exclusions]
1       2       3       4
2       1       2       3
3       1       2       4
4       1       2       3

[dummies3]
; Dummy from            funct   a       b
4	1	2	3	1	0.13458335      0.13458335

[system]
iceIh

[molecules]
S   1
water 439 

EOF 
    
    #lambda=`echo "${k_spring}"|awk '{print $1/2}'` #U_spring=lambda*(r-r0)^2; k=2*lambda

    for i in `seq 0 1 $length`
    do
    mkdir ${Temperature[${i}]}_${Pressure[${i}]}
    mkdir ${Temperature[${i}]}_${Pressure[${i}]}/Einstein
    #minimization
    mkdir ${Temperature[${i}]}_${Pressure[${i}]}/em
    nmol=`grep OW ${Temperature[${i}]}_${Pressure[${i}]}/eq/eq.gro|wc -l`
    gmx grompp -f ./mdp/em.mdp -c ${Temperature[${i}]}_${Pressure[${i}]}/eq/eq.gro -p ${Temperature[${i}]}_${Pressure[${i}]}/eq/topo.top -o ${Temperature[${i}]}_${Pressure[${i}]}/em/em
    gmx mdrun -v -deffnm ${Temperature[${i}]}_${Pressure[${i}]}/em/em
    echo "0"|gmx trjconv -f ${Temperature[${i}]}_${Pressure[${i}]}/em/em.trr -s ${Temperature[${i}]}_${Pressure[${i}]}/em/em.tpr -dump -1 -o ${Temperature[${i}]}_${Pressure[${i}]}/em/em.g96

    #mdp
    cat mdp/einstein.mdp >${Temperature[${i}]}_${Pressure[${i}]}/Einstein/einstein.mdp
    sed -i s"@ref_t.*@ref_t=${Temperature[${i}]}@" ${Temperature[${i}]}_${Pressure[${i}]}/Einstein/einstein.mdp
    sed -i s"@gen_temp.*@gen_temp=${Temperature[${i}]}@" ${Temperature[${i}]}_${Pressure[${i}]}/Einstein/einstein.mdp
    sed -i s"@coulombtype.*@coulombtype=cut-off@" ${Temperature[${i}]}_${Pressure[${i}]}/Einstein/einstein.mdp
    #topo
    k_spring=`echo "${lambda0} ${Temperature[${i}]}"|awk '{print $1*$2*2}'` #U_spring=0.5k_spring(r-r0)^2
    sed '$d' ./src/ideal.template >${Temperature[${i}]}_${Pressure[${i}]}/Einstein/ideal.top
    sed -i s"@spring_const.@${k_spring}@"g ${Temperature[${i}]}_${Pressure[${i}]}/Einstein/ideal.top
    echo "water $nmol"|awk '{print $1, $2-1}' >>${Temperature[${i}]}_${Pressure[${i}]}/Einstein/ideal.top
    #index
    echo -e "r 1 & a OW\n name 3 SWF1 \n q"|gmx make_ndx -f ${Temperature[${i}]}_${Pressure[${i}]}/em/em.gro -o ${Temperature[${i}]}_${Pressure[${i}]}/index.ndx

    gmx grompp -f ${Temperature[${i}]}_${Pressure[${i}]}/Einstein/einstein.mdp -c ${Temperature[${i}]}_${Pressure[${i}]}/em/em.g96 -p ${Temperature[${i}]}_${Pressure[${i}]}/Einstein/ideal.top -o ${Temperature[${i}]}_${Pressure[${i}]}/Einstein/eq -n ${Temperature[${i}]}_${Pressure[${i}]}/index.ndx -r ${Temperature[${i}]}_${Pressure[${i}]}/em/em.g96 -maxwarn 1




cat <<EOF >eq.sh 
#!/bin/bash
#SBATCH -n 16
#SBATCH -J ${Temperature[${i}]}_${Pressure[${i}]}_einstein
#SBATCH -N 1
#SBATCH -p NVIDIAGeForceRTX3090
#SBATCH --gres=gpu:1
#SBATCH -e error.err1
#SBATCH -o output.out1

module purge
module load compiler/gcc/7.3.1
module load compiler/intel/2021.3.0
module load mpi/intelmpi/2021.3.0
module load mathlib/fftw/intelmpi/3.3.9_single
module load apps/gromacs/intelmpi/2018.8-30
gmx mdrun -pin on -deffnm ./${Temperature[${i}]}_${Pressure[${i}]}/Einstein/eq -v #-gpu_id 0 -pme gpu -nb gpu
EOF

sbatch eq.sh
    done
fi


#Decoupl
if [ "${Decoupl}" = "yes" ]; then
    
cat <<EOF >cpt_lambda.py
from math import *
import sys, random,os,shutil
from scipy.integrate import quad
from scipy.special import roots_legendre
import numpy as np
T=float(sys.argv[1])
k_max=T*5000*2*8.314/10#831400#
print(f'k_max={k_max}')
num =int(sys.argv[2])
c = np.exp(3.5)#(9.30685)#(3.5)#(9.30685)##(8.6225)##
k=[]
points,weights=roots_legendre(num)

#print(points)
for i in range(len(points)):
    points[i]=(points[i]+1)/2
#print(points)
x_i=[]
E_i=[]
for i in range(num):
    x_i.append(log(c)+(log(k_max)-log(c))*points[i])
    E_i.append(exp(x_i[i])-c)
#print(x_i)
#print(E_i)
for i in range(num):
   k.append(exp(log(c)+points[i]*(log(k_max+c)-log(c)) )-c)
print('\nspring constant:\n')
print(k)
print('\nweights:\n')
print(weights[:])


with open('lambda.txt','w') as f:
    print( np.round(k_max,5), np.round(c,5), np.round(c/2*(np.log(k_max+c)-np.log(c)),5), file=f)
    for i in range(len(k)):
        print( np.round(k[i],5), np.round(weights[i],5), np.round(weights[i]/2*(np.log(k_max+c)-np.log(c)),5), file=f)

EOF


cat <<EOF >src/decoupl.template
[ defaults ]
; nbfunc    comb-rule   gen-pairs   fudgeLJ fudgeQQ
  1     2       no      1.0 1.0

[atomtypes]
;name     mass      charge   ptype    sigma        epsilon
MW     0             0.000       D   0.0           0.0
OW     15.99940      0.000       A   3.16680e-01  8.82164e-01
HW     1.00800       0.000       A   0.00000E+00   0.00000E+00
SW     15.99940      0.000       A   3.16680e-01  8.82164e-01

[ moleculetype ]
; molname       nrexcl
S                1

[atoms]
; nr type resnr residu atom cgnr charge
1     SW   1     S  SWF1  1     0         15.994
2     HW   1     S  HWF2  1     0.5897    1.008
3     HW   1     S  HWF3  1     0.5897    1.008
4     MW   1     S  MWF4  1    -1.1794    0.0

[ position_restraints ]
; ai  funct  fcx    fcy    fcz
   2    1    spring_const. spring_const. spring_const. ; restrains to a point
   3    1    spring_const. spring_const. spring_const. ; restrains to a point

[constraints]
;i j funct doh  dhh
1       2       1       0.09572
1       3       1       0.09572
2       3       1       0.15139

[exclusions]
1       2       3       4
2       1       2       3
3       1       2       4
4       1       2       3

[dummies3]
; Dummy from                    funct   a               b
4	1	2	3	1	0.13458335      0.13458335

[moleculetype]
; name nrexcl
water  1

[atoms]
; nr type resnr residu atom cgnr charge
1     OW   1     water  OW1  1     0         15.994 
2     HW   1     water  HW2  1     0.5897    1.008
3     HW   1     water  HW3  1     0.5897    1.008
4     MW   1     water  MW4  1    -1.1794    0.0

[ position_restraints ]
; ai  funct  fcx    fcy    fcz
   1    1    spring_const. spring_const. spring_const. ; restrains to a point
   2    1    spring_const. spring_const. spring_const. ; restrains to a point
   3    1    spring_const. spring_const. spring_const. ; restrains to a point

[constraints]
;i j funct doh  dhh
1       2       1       0.09572
1       3       1       0.09572
2       3       1       0.15139

[exclusions]
1       2       3       4
2       1       2       3
3       1       2       4
4       1       2       3

[dummies3]
; Dummy from            funct   a       b
4	1	2	3	1	0.13458335      0.13458335

[system]
iceIh

[molecules]
S   1
water 639 

EOF
    
    #lambda=`echo "${k_spring}"|awk '{print $1/2}'` #U_spring=lambda*(r-r0)^2; k=2*lambda

    for i in `seq 0 1 $length`
    do
    mkdir ${Temperature[${i}]}_${Pressure[${i}]}
    mkdir ${Temperature[${i}]}_${Pressure[${i}]}/Decoupl
    nmol=`grep OW ${Temperature[${i}]}_${Pressure[${i}]}/eq/eq.gro|wc -l`

    #mdp
    cat mdp/einstein.mdp >${Temperature[${i}]}_${Pressure[${i}]}/Decoupl/einstein.mdp
    sed -i s"@ref_t.*@ref_t=${Temperature[${i}]}@" ${Temperature[${i}]}_${Pressure[${i}]}/Decoupl/einstein.mdp
    sed -i s"@gen_temp.*@gen_temp=${Temperature[${i}]}@" ${Temperature[${i}]}_${Pressure[${i}]}/Decoupl/einstein.mdp
    #topo
        #for f in `sed -n "2,$"p src/landa_weight.dat|awk '{print $2}'`
        #do
        #mkdir ${Temperature[${i}]}_${Pressure[${i}]}/Decoupl/${f}
        #k_spring=`echo "${lambda0} ${Temperature[${i}]} ${f}"|awk '{print $1*$2*2*$3}'` #scaled k
        #sed '$d' ./src/decoupl.template >${Temperature[${i}]}_${Pressure[${i}]}/Decoupl/${f}/decoupl.top
        #sed -i s"@spring_const.@${k_spring}@"g ${Temperature[${i}]}_${Pressure[${i}]}/Decoupl/${f}/decoupl.top
        #echo "water $nmol"|awk '{print $1, $2-1}' >>${Temperature[${i}]}_${Pressure[${i}]}/Decoupl/${f}/decoupl.top
        
        python3 cpt_lambda.py ${Temperature[${i}]} ${num_points} 
        mv lambda.txt ./${Temperature[${i}]}_${Pressure[${i}]}/lambda.txt 
        for f in `sed -n "2,$"p ${Temperature[${i}]}_${Pressure[${i}]}/lambda.txt|awk '{print $1}'`
#        for f in `sed -n "2"p src/lambda.txt|awk '{print $1}'`
        do
        mkdir ${Temperature[${i}]}_${Pressure[${i}]}/Decoupl/${f}
        k_spring=${f} #scaled k
        sed '$d' ./src/decoupl.template >${Temperature[${i}]}_${Pressure[${i}]}/Decoupl/${f}/decoupl.top
        sed -i s"@spring_const.@${k_spring}@"g ${Temperature[${i}]}_${Pressure[${i}]}/Decoupl/${f}/decoupl.top
        echo "water $nmol"|awk '{print $1, $2-1}' >>${Temperature[${i}]}_${Pressure[${i}]}/Decoupl/${f}/decoupl.top

        gmx grompp -f ${Temperature[${i}]}_${Pressure[${i}]}/Decoupl/einstein.mdp -c ${Temperature[${i}]}_${Pressure[${i}]}/em/em.g96 -p ${Temperature[${i}]}_${Pressure[${i}]}/Decoupl/${f}/decoupl.top -o ${Temperature[${i}]}_${Pressure[${i}]}/Decoupl/${f}/eq -n ${Temperature[${i}]}_${Pressure[${i}]}/index.ndx -r ${Temperature[${i}]}_${Pressure[${i}]}/em/em.g96 -maxwarn 1



cat <<EOF >eq.sh 
#!/bin/bash
#SBATCH -n 16
#SBATCH -J decoupl_${Temperature[${i}]}_${Pressure[${i}]}_${f}
#SBATCH -N 1
#SBATCH -p NVIDIAGeForceRTX3090
#SBATCH --gres=gpu:1
#SBATCH -e error.err1
#SBATCH -o output.out1

module purge
module load compiler/gcc/7.3.1
module load compiler/intel/2021.3.0
module load mpi/intelmpi/2021.3.0
module load mathlib/fftw/intelmpi/3.3.9_single
module load apps/gromacs/intelmpi/2018.8-30
gmx mdrun -ntmpi 1 -ntomp 16 -pin on -deffnm ./${Temperature[${i}]}_${Pressure[${i}]}/Decoupl/${f}/eq -v -gpu_id 0 -pme gpu -nb gpu
EOF

sbatch eq.sh
        done
    done
fi

if [ "${Analyze_Einstein}" = "yes" ]; then
    for i in `seq 0 1 $length`
    do
    echo "T=${Temperature[${i}]} K, p=${Pressure[${i}]} bar" >> ${Temperature[${i}]}_${Pressure[${i}]}/result.dat
    mkdir ${Temperature[${i}]}_${Pressure[${i}]}/Einstein/rerun
    nmol=`grep OW ${Temperature[${i}]}_${Pressure[${i}]}/eq/eq.gro|wc -l`
    gmx mdrun -v -deffnm ${Temperature[${i}]}_${Pressure[${i}]}/Einstein/rerun/rerun -s ${Temperature[${i}]}_${Pressure[${i}]}/eq/eq.tpr -rerun ${Temperature[${i}]}_${Pressure[${i}]}/Einstein/eq.xtc
    echo -e "Potential \n 0"|gmx energy -f ${Temperature[${i}]}_${Pressure[${i}]}/Einstein/rerun/rerun.edr -b 1000 -o ${Temperature[${i}]}_${Pressure[${i}]}/Einstein/rerun/rerun.xvg
    tail ${Temperature[${i}]}_${Pressure[${i}]}/Einstein/rerun/rerun.xvg -n 300 |awk '{print $2*1000/'${nmol}'/8.314/'${Temperature[${i}]}'}'>${Temperature[${i}]}_${Pressure[${i}]}/Einstein/Usol.dat
    
cat <<EOF >src/A1.py
import numpy as np
import sys

path=sys.argv[1]
Usol=np.loadtxt('{}/Usol.dat'.format(path))
Ulatt=Usol[0]

A1=Ulatt-np.log(np.mean(np.exp(Ulatt-Usol)))
print('A1(NkT):', A1, ",  Usol(NkT):", np.mean(Usol))
EOF
    python3 ./src/A1.py "${Temperature[${i}]}_${Pressure[${i}]}/Einstein" >>${Temperature[${i}]}_${Pressure[${i}]}/result.dat

    done

fi

if [ "${Analyze_Decoupl}" = "yes" ]; then

cat <<EOF >src/A2.py 
import numpy as np
import sys

T=sys.argv[1]
P=sys.argv[2]
U_spr=np.loadtxt(f'{T}_{P}/U_spring.dat')
k=np.loadtxt(f'{T}_{P}/lambda.txt')

result=np.sum(U_spr/k[1:,0]*(k[1:,0]+np.exp(3.5))*k[1:,2])*1000/8.314/float(T)

print(f'A2(NkT): {result}')

EOF
    for i in `seq 0 1 $length`
    do
    nmol=`grep OW ${Temperature[${i}]}_${Pressure[${i}]}/eq/eq.gro|wc -l`
    #echo "#lambda MSD" >${Temperature[${i}]}_${Pressure[${i}]}/Decoupl/landa_msd.dat
    echo "#U_spring kJ/mol" >${Temperature[${i}]}_${Pressure[${i}]}/U_spring.dat
        #for f in `sed -n "2,$"p src/landa_weight.dat|awk '{print $2}'`
        for f in `sed -n "2,$"p ${Temperature[${i}]}_${Pressure[${i}]}/lambda.txt|awk '{print $1}'`
        do
        k_spring=`echo "${lambda0} ${Temperature[${i}]} ${f}"|awk '{print $1*$2*2*$3}'` #scaled k
        lambda=`echo "${k_spring}"|awk '{print $1/2}'` #U_spring=lambda*(r-r0)^2; k=2*lambda
        U_spring2=`echo -e "Position-Rest. \n 0"|gmx energy -f ${Temperature[${i}]}_${Pressure[${i}]}/Decoupl/${f}/eq.edr -b 1000 -o ${Temperature[${i}]}_${Pressure[${i}]}//Decoupl/${f}/posre.xvg |grep Position|awk '{print $3/'${nmol}'}'` 
        echo "${U_spring2}" >>${Temperature[${i}]}_${Pressure[${i}]}/U_spring.dat
#        U_spring=`echo -e "Position-Rest. \n 0"|gmx energy -f ${Temperature[${i}]}_${Pressure[${i}]}/Decoupl/${f}/eq.edr -b 1000 -o ${Temperature[${i}]}_${Pressure[${i}]}//Decoupl/${f}/posre.xvg |grep Position|awk '{print $3*1000/'${nmol}'/8.314/'${Temperature[${i}]}'}'`
#        echo "${f} ${U_spring} ${lambda}"|awk '{print $1, $2/$3}' >>${Temperature[${i}]}_${Pressure[${i}]}/Decoupl/landa_msd.dat
#        T_sim=`echo -e "Temperature \n 0"|gmx energy -f ${Temperature[${i}]}_${Pressure[${i}]}/Decoupl/${f}/eq.edr -b 1000 -o ${Temperature[${i}]}_${Pressure[${i}]}//Decoupl/${f}/temp.xvg|grep Temperature|awk '{print $2}'`
#        echo "${f} ${T_sim}" >>${Temperature[${i}]}_${Pressure[${i}]}/Decoupl/landa_temp.dat
        done

#    deltaA2=`paste ${Temperature[${i}]}_${Pressure[${i}]}/U_spring.dat ${Temperature[${i}]}_${Pressure[${i}]}/lambda.txt |sed -n "2,$"p |awk '{print ($1/$2)*($1+exp(3.5))*$4*1000/8.314/'${Temperature[${i}]}' }'|awk '{total+=$1} END {print -1*total}'`

#    echo "A2(NkT): ${deltaA2}" >>${Temperature[${i}]}_${Pressure[${i}]}/result.dat
    python3 ./src/A2.py ${Temperature[${i}]} ${Pressure[${i}]} >>${Temperature[${i}]}_${Pressure[${i}]}/result.dat
    echo "A_symm(NkT): -0.69" >>${Temperature[${i}]}_${Pressure[${i}]}/result.dat
    echo "A_pauling(NkT): -0.405" >>${Temperature[${i}]}_${Pressure[${i}]}/result.dat
    pV=`echo -e "pV \n 0"|gmx energy -f ${Temperature[${i}]}_${Pressure[${i}]}/eq/eq.edr -o ${Temperature[${i}]}_${Pressure[${i}]}/eq/pV.xvg|grep pV|awk '{print $2*1000/8.314/'${Temperature[${i}]}'/'${nmol}'}'`
    echo "pV(NkT): ${pV}" >> ${Temperature[${i}]}_${Pressure[${i}]}/result.dat

    done

fi
