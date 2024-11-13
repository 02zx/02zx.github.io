#2024/11/13/24:00-ver
#settings

Pressure=("none" "none" "none" "none" "none" "none")    #bar
Temperature=("200" "250" "300" "350" "400" "450")  #K 

N_site=4 #3 for 3 site model, 4 for 4 site model.
water_model=tip4p #spce/tip4p/tip4p2005/tip4pice...

create_system=yes
N=360
density=1.05 #g/cm^3
box_size=1.9
conf=conf
top=topo
state=liquid #solid or liquid

Minimization=yes
Equilibration=yes

#--Hamiltonian Integration---
#
Integration=no
num_points=14 #number of point for integration
Analyze=no #compute integral
#----------------------------

#--Thermodynamic Integration---
#
TI=yes
#------------------------------ 




#---------------------------------------------------------------------
let length=${#Pressure[@]}-1
ensemble=npt
mkdir src

cat<<EOF > topo.top
Include forcefield parameters
#include "amber99sb.ff/forcefield.itp"
#include "./src/water_model.itp"

[ system ]
; Name
ICE in water

[ molecules ]
; Compound        #mols
EOF

if [ "${water_model}" = "spce" ]; then

    charge=0.4238
    sigma=3.1656 #Angstrom
    epsilon=78.20 #K
    cat<<EOF >./src/spce.template
[ moleculetype ]
; molname	nrexcl
SOL		2

[ atoms ]
; id  at type     res nr  res name  at name  cg nr  charge    mass
  1   OW_spc      1       SOL       OW       1      QO   15.99940
  2   HW_spc      1       SOL       HW1      1      QH    1.00800
  3   HW_spc      1       SOL       HW2      1      QH    1.00800

#ifndef FLEXIBLE

[ settles ]
; OW	funct	doh	dhh
1       1       0.1     0.16330

[ exclusions ]
1	2	3
2	1	3
3	1	2

#else

[ bonds ]
; i     j       funct   length  force.c.
1       2       1       0.1     345000  0.1     345000
1       3       1       0.1     345000  0.1     345000

[ angles ]
; i     j       k       funct   angle   force.c.
2       1       3       1       109.47  383     109.47  383

#endif
EOF

elif [ "${water_model}" = "tip4p" ]; then

    charge=0.52
    sigma=3.154
    epsilon=78.0
    cat<<EOF >./src/tip4p.template
[ moleculetype ]
; molname	nrexcl
SOL		2

[ atoms ]
; id  at type     res nr  res name  at name  cg nr  charge    mass
  1   OW_tip4p    1       SOL       OW       1       0        16.00000
  2   HW_tip4p    1       SOL       HW1      1       QH      1.00800
  3   HW_tip4p    1       SOL       HW2      1       QH      1.00800
  4   MW          1       SOL       MW       1       QO      0.00000

#ifndef FLEXIBLE

[ settles ]
; i	funct	doh	dhh
1	1	0.09572	0.15139

#else

[ bonds ]
; i     j       funct   length  force.c.
1       2       1       0.09572 502416.0 0.09572        502416.0 
1       3       1       0.09572 502416.0 0.09572        502416.0 
        
[ angles ]
; i     j       k       funct   angle   force.c.
2       1       3       1       104.52  628.02  104.52  628.02  

#endif


[ virtual_sites3 ]
; Vsite from                    funct   a               b
4       1       2       3       1       0.128012065     0.128012065


[ exclusions ]
1	2	3	4
2	1	3	4
3	1	2	4
4	1	2	3


; The position of the virtual site is computed as follows:
;
;		O
;  	      
;	    	V
;	  
;	H		H
;
; const = distance (OV) / [ cos (angle(VOH)) 	* distance (OH) ]
;	  0.015 nm	/ [ cos (52.26 deg)	* 0.09572 nm	]
;
; Vsite pos x4 = x1 + a*(x2-x1) + b*(x3-x1)
EOF

elif [ "${water_model}" = "tip4p2005" ]; then

    charge=0.5564
    sigma=3.1589
    epsilon=93.2
    cat<<EOF >./src/tip4p2005.template
[ atomtypes ]
; tip4p_2005

  HW_2005     1       1.008   0.0000  A   0.00000e+00  0.00000e+00

  OW_2005     8       16.00   0.0000  A   3.15890e-01  7.74898e-01

  MW          0      0.0000   0.0000  D   0.00000e+00  0.00000e+00


[ moleculetype ]
; molname	nrexcl
SOL		2

[ atoms ]
; id  at type     res nr  res name  at name  cg nr  charge    mass
  1   OW_2005    1       SOL       OW       1       0        16.00000
  2   HW_2005    1       SOL       HW1      2       QH      1.00800
  3   HW_2005    1       SOL       HW2      3       QH    1.00800
  4   MW          1       SOL       MW      4       QO     0.00000

#ifndef FLEXIBLE

[ settles ]
; i	funct	doh	dhh
1	1	0.09572	0.15139

#else

[ bonds ]
; i     j       funct   length  force.c.
1       2       1       0.09572 502416.0 0.09572        502416.0 
1       3       1       0.09572 502416.0 0.09572        502416.0 
        
[ angles ]
; i     j       k       funct   angle   force.c.
2       1       3       1       104.52  628.02  104.52  628.02  

#endif


[ virtual_sites3 ]
; Vsite from                    funct   a               b
4       1       2       3       1       0.13193828     0.13193828


[ exclusions ]
1	2	3	4
2	1	3	4
3	1	2	4
4	1	2	3


; The position of the virtual site is computed as follows:
;
;		O
;  	      
;	    	V
;	  
;	H		H
;
; const = distance (OV) / [ cos (angle(VOH)) 	* distance (OH) ]
;	  0.015 nm	/ [ cos (52.26 deg)	* 0.09572 nm	]
;
; Vsite pos x4 = x1 + a*(x2-x1) + b*(x3-x1)
EOF

elif [ "${water_model}" = "tip4pice" ]; then

    charge=0.5897
    sigma=3.1668
    epsilon=106.1
    cat<<EOF >./src/tip4pice.template
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
2       H_ice        1       SOL      HW1     1       QH
3       H_ice        1       SOL      HW2     1       QH
4       M            1       SOL      MW      1       QO

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
EOF

fi 

qh=${charge} 
qo=`echo "${qh}"|awk '{print -1*$qh*2}'`
cat ./src/${water_model}.template >./src/${water_model}.itp
sed -i s"@QO@$qo@" ./src/${water_model}.itp
sed -i s"@QH@$qh@"g ./src/${water_model}.itp
sed -i s"@water_model.itp@${water_model}.itp@" ${top}.top 
if [ "`tail -n 1 topo.top |awk '{print $1}'`" = "SOL" ]; then

    sed -i s"@^SOL.*@SOL ${N}@" ${top}.top

else

    echo "SOL $N" >> ${top}.top

fi

#generate simulation box
if [ "${create_system}" = "yes" ]; then
    if [ "${N_site}" = "3" ]; then

        cat<<EOF > water.pdb
TITLE     Generated by gmx solvate
REMARK    THIS IS A SIMULATION BOX
MODEL        1
ATOM      1  OW  SOL     1       2.300   6.280   1.130  1.00  0.00
ATOM      2  HW1 SOL     1       1.370   6.260   1.500  1.00  0.00
ATOM      3  HW2 SOL     1       2.310   5.890   0.210  1.00  0.00
EOF


    fi
    
    if [ "${N_site}" = "4" ]; then
        cat<<EOF > water.pdb
TITLE     Generated by gmx solvate
REMARK    THIS IS A SIMULATION BOX
MODEL        1
ATOM      1  OW  SOL     1      17.360   8.390   2.570  1.00  0.00
ATOM      2  HW1 SOL     1      17.770   7.810   3.220  1.00  0.00
ATOM      3  HW2 SOL     1      16.430   8.310   2.740  1.00  0.00
ATOM      4  MW  SOL     1      17.300   8.310   2.670  1.00  0.00
EOF

    fi

    if [ "${density}" != "no" ]; then 
    
        box_size=`echo "${N} ${density}"|awk '{printf "%.3f" ,($1*18/602/$2)^(1/3)}'`    

    fi
    gmx insert-molecules -box ${box_size} -ci water.pdb -nmol ${N} -try 2000 -o ${conf}.pdb 
fi
    
#Minimization
if [ "${Minimization}" = "yes" ]; then
    mkdir mdp
    mkdir em
    cat<<EOF >mdp/em.mdp
;define = -DFLEXIBLE
integrator = steep
nsteps = 10000
emtol  = 10.0
emstep = 0.01
;
nstxout   = 100
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
    gmx grompp -f mdp/em.mdp -c ${conf}.pdb -p ${top} -o em/em
    gmx mdrun -v -deffnm em/em
    rm step*

fi

#Equilibration
if [ "${Equilibration}" = "yes" ]; then

    cat<<EOF >>simulation.log
#---Equilibration----
Thermodynamic Integration: ${TI}
Ensemble:$ensemble"
Density: $density g/cm^3"
Box Size: $box_size nm"
Number of Water: $N"
Water model: ${water_model}"
EOF
    cat<<'EOF' >eq.sh
#!/bin/bash
#SBATCH -n 16
#SBATCH -J equilibration
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
module load apps/gromacs/2023.4-30

T=
P=
gmx mdrun -ntmpi 1 -ntomp 16 -pin on -deffnm ./${T}_${P}/eq/eq -v -gpu_id 0 -pme gpu -nb gpu
EOF


    
    mkdir mdp
    cat<<EOF >mdp/eq.mdp
;VARIOUS PREPROCESSING OPTIONS
title                    = equilibration

; RUN CONTROL PARAMETERS
integrator               = md
; Start time and timestep in ps
tinit                    = 0
dt                       = 0.002
nsteps                   = 1000000
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
tau_t                    = 2
ref_t=443
tc-grps                  = system
EOF
    cp mdp/eq.mdp mdp/md.mdp
    sed -i s"@^dt.*@dt=0.001@" mdp/md.mdp 
    sed -i s"@^nsteps.*@nsteps=5000000@" mdp/md.mdp 
    cat<<EOF >> mdp/md.mdp  
; GENERATE VELOCITIES FOR STARTUP RUN
 gen_vel                  = yes
 gen_temp=443
 gen_seed                 = 50934891

; OPTIONS FOR BONDS
constraints              = hbonds ;all-angles
; Type of constraint algorithm
constraint-algorithm     = lincs
lincs-iter               =  4
lincs-order              =  6
EOF

    if [ "${density}" != "no" ]; then

        ensemble=nvt

    else 

        cat<<EOF >> mdp/eq.mdp
Pcoupl              =  c-rescale
Pcoupltype          =  isotropic
tau_p               =  1       
compressibility     =  4.5e-5 
ref_p=4010
EOF

    fi

    cat<<EOF >> mdp/eq.mdp  
; GENERATE VELOCITIES FOR STARTUP RUN
 gen_vel                  = yes
 gen_temp=443
 gen_seed                 = 50934891

; OPTIONS FOR BONDS
constraints              = hbonds ;all-angles
; Type of constraint algorithm
constraint-algorithm     = lincs
lincs-iter               =  4
lincs-order              =  6
EOF
    if [ "${TI}" = "yes" ]; then 
        sed -i s"@^nsteps.*@nsteps=2500000@" mdp/eq.mdp 
    fi
    
    for i in `seq 0 1 $length`
    do

    if [ "${density}" != "no" ]; then 
        Pressure[${i}]='none'
        ensemble=nvt
    fi
       mkdir ${Temperature[${i}]}_${Pressure[${i}]}
       mkdir ${Temperature[${i}]}_${Pressure[${i}]}/eq
       sed -i s"@ref_t.*@ref_t=${Temperature[${i}]}@" mdp/eq.mdp
       sed -i s"@gen_temp.*@gen_temp=${Temperature[${i}]}@" mdp/eq.mdp

    if [ "${state}" = "solid" ]; then 
        sed -i s"@^Pcouple.*@Pcouple = Berendsen@" mdp/eq.mdp 
        sed -i s"@^Pcoupltype.*@Pcoupltype = anisotropic@" mdp/eq.mdp 
        sed -i s"@^compressibility.*@compressibility = 4.5e-5 4.5e-5 4.5e-5 0 0 0 @" mdp/eq.mdp 
        sed -i s"@^ref_p.*@ref_p =${Pressure[${i}]} ${Pressure[${i}]} ${Pressure[${i}]} 0 0 0  @" mdp/eq.mdp 
    else
       sed -i s"@ref_p.*@ref_p=${Pressure[${i}]}@" mdp/eq.mdp
    fi

       gmx grompp -f mdp/eq.mdp -c em/em.gro -p ${top}.top -o ${Temperature[${i}]}_${Pressure[${i}]}/eq/eq  
   
       sed -i s"@^T=.*@T=${Temperature[${i}]}@" eq.sh
       sed -i s"@^P=.*@P=${Pressure[${i}]}@" eq.sh
       sbatch eq.sh
       echo "Temperature:${Temperature[${i}]} K, Pressure: ${Pressure[${i}]} bar" >>simulation.log
    done
fi

#Integration
if [ "${Integration}" = "yes" ]; then
    cat<<EOF >>simulation.log
#---Integration----
Ensemble:$ensemble"
Density: $density g/cm^3"
Box Size: $box_size nm"
Number of Water: $N"
Water model: ${water_model}"
EOF


    mkdir src mdp
    cat<<EOF >src/gauss_legendre.py
import numpy as np
import sys

num=int(sys.argv[1])
nodes, weights = np.polynomial.legendre.leggauss(num) 
nodes_transformed = 0.5 * (nodes + 1)
weights_transformed = weights * 0.5
print('#sqrt(lambda), lambda, weight')
for i in range(num):
    print(np.round(np.sqrt(nodes_transformed[i]),5),np.round(nodes_transformed[i],5),np.round(weights_transformed[i],5))

#integral = np.sum(weights_transformed * f(nodes_transformed))

#print("The integral from 0 to 1 is:", integral)
EOF

    python3 src/gauss_legendre.py $num_points >src/landa_weight.dat
    gmx grompp -f mdp/md.mdp -c em/em.gro -p ${top}.top -o full.tpr
    for i in `seq 0 1 $length`
    do
        if [ "${density}" != "no" ]; then 
            Pressure[${i}]='none'
            ensemble=nvt
        fi

        
        for f in `sed -n "2,$"p src/landa_weight.dat|awk '{print $1}'` #loop for sqrt(lambda)
        do
            mkdir ${Temperature[${i}]}_${Pressure[${i}]}/$f
            qh=`echo "${charge} $f"|awk '{printf "%.5f", $1*$2}'` 
            qo=`echo "${qh}"|awk '{print -1*$qh*2}'`
            
            cat ./src/${water_model}.template >${Temperature[${i}]}_${Pressure[${i}]}/$f/${water_model}.itp
            sed -i s"@QO@$qo@" ${Temperature[${i}]}_${Pressure[${i}]}/$f/${water_model}.itp
            sed -i s"@QH@$qh@"g ${Temperature[${i}]}_${Pressure[${i}]}/$f/${water_model}.itp
            cat ${top}.top|sed s"@src/${water_model}.itp@${water_model}.itp@" >${Temperature[${i}]}_${Pressure[${i}]}/$f/${top}.top
            
            #minimization
            mkdir ${Temperature[${i}]}_${Pressure[${i}]}/$f/em
            gmx grompp -f mdp/em.mdp -c ${Temperature[${i}]}_${Pressure[${i}]}/eq/eq.gro -p ${Temperature[${i}]}_${Pressure[${i}]}/$f/${top}.top -o ${Temperature[${i}]}_${Pressure[${i}]}/$f/em/em
            gmx mdrun -v -deffnm ${Temperature[${i}]}_${Pressure[${i}]}/$f/em/em
            rm step*

            #nvt simulation
            mkdir ${Temperature[${i}]}_${Pressure[${i}]}/$f/md
            sed -i s"@ref_t.*@ref_t=${Temperature[${i}]}@" mdp/md.mdp
            gmx grompp -f mdp/md.mdp -c ${Temperature[${i}]}_${Pressure[${i}]}/$f/em/em.gro -p ${Temperature[${i}]}_${Pressure[${i}]}/$f/${top}.top -o ${Temperature[${i}]}_${Pressure[${i}]}/$f/md/md
            cat<<'EOF' >md.sh
#!/bin/bash
#SBATCH -n 16
#SBATCH -J nvt
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
module load apps/gromacs/2023.4-30

T=
P=
factor=
gmx mdrun -ntmpi 1 -ntomp 16 -pin on -deffnm ./${T}_${P}/${factor}/md/md -v -gpu_id 0 -pme gpu -nb gpu
mkdir ${T}_${P}/${factor}/rerun

gmx mdrun -v -s full.tpr -deffnm ${T}_${P}/${factor}/rerun/rerun -rerun ./${T}_${P}/${factor}/md/md.xtc 
EOF
            sed -i s"@factor=.*@factor=${f}@" md.sh
            sed -i s"@^T=.*@T=${Temperature[${i}]}@" md.sh
            sed -i s"@^P=.*@P=${Pressure[${i}]}@" md.sh
            sbatch md.sh
            echo "Temperature:${Temperature[${i}]} K, Pressure: ${Pressure[${i}]} bar" >>simulation.log

        done
    done
fi


#Integration
#"${state} ="Solid" &&  "${TI}" = "yes" 
#if [ "${TI}" = "yes" ]; then



#fi



#Analyze

if [ "${Analyze}" = "yes" ]; then #TI=no
cat<<'EOF' >MBWR.py
import numpy as np
import sys
    
    

def MBWR(T,rho):
    #T and rho are in reduced units.

    #Johnson J K, Zollweg J A and Gubbins K E 1993 Mol.Phys. 78 591
    #See Tab.3 for validation
    
    #Tab.10
    x1=0.8623085097507421
    x2=2.976218765822098
    x3=-8.402230115796038
    x4=0.1054136629203555
    x5=-0.8564583828174598
    x6=1.582759470107601
    x7=0.7639421948305453
    x8=1.753173414312048
    x9=2.798291772190376e+3
    x10=-4.8394220260857657e-2
    x11=0.9963265197721935
    x12=-3.698000291272493e+1
    x13=2.084012299434647e+1
    x14=8.305402124717285e+1
    x15=-9.574799715203068e+2
    x16=-1.477746229234994e+2
    x17=6.398607852471505e+1
    x18=1.603993673294834e+1
    x19=6.805916615864377e+1
    x20=-2.791293578795945e+3
    x21=-6.245128304568454
    x22=-8.116836104958410e+3
    x23=1.488735559561229e+1
    x24=-1.059346754655084e+4
    x25=-1.131607632802822e+2
    x26=-8.867771540418822e+3
    x27=-3.986982844450543e+1
    x28=-4.689270299917261e+3
    x29=2.593535277438717e+2
    x30=-2.694523589434903e+3
    x31=-7.218487631550215e+2
    x32=1.721802063863269e+2
    
    #Tab. 5
    a1=x1*T+x2*np.sqrt(T)+x3+x4/T+x5/T**2
    a2=x6*T+x7+x8/T+x9/T**2
    a3=x10*T+x11+x12/T
    a4=x13
    a5=x14/T+x15/T**2
    a6=x16/T
    a7=x17/T+x18/T**2
    a8=x19/T**2
    
    #Tab.6
    b1=x20/T**2+x21/T**3
    b2=x22/T**2+x23/T**4
    b3=x24/T**2+x25/T**3
    b4=x26/T**2+x27/T**4
    b5=x28/T**2+x29/T**3
    b6=x30/T**2+x31/T**3+x32/T**4
    
    #Tab.7
    gamma=3
    F=np.exp(-gamma*rho**2)
    G1=(1-F)/(2*gamma)
    G2=-(F*rho**2-2*G1)/(2*gamma)
    G3=-(F*rho**4-4*G2)/(2*gamma)
    G4=-(F*rho**6-6*G3)/(2*gamma)
    G5=-(F*rho**8-8*G4)/(2*gamma)
    G6=-(F*rho**10-10*G5)/(2*gamma)
    
    #Tab.8
    c1=x2*np.sqrt(T)/2+x3+2*x4/T+3*x5/T**2
    c2=x7+2*x8/T+3*x9/T**2
    c3=x11+2*x12/T
    c4=x13
    c5=2*x14/T+3*x15/T**2
    c6=2*x16/T
    c7=2*x17/T+3*x18/T**2
    c8=3*x19/T**2
    
    #Tab.9
    d1=3*x20/T**2+4*x21/T**3
    d2=3*x22/T**2+5*x23/T**4
    d3=3*x24/T**2+4*x25/T**3
    d4=3*x26/T**2+5*x27/T**4
    d5=3*x28/T**2+4*x29/T**3
    d6=3*x30/T**2+4*x31/T**3+5*x32/T**4
    


        

    #Ar(NVT)=A(NVT)-Aid(NVT); Ar_ = Ar/N/epsilon
    #eq.5
    Ar_=(a1*rho**1)/1 \
      + (a2*rho**2)/2 \
      + (a3*rho**3)/3 \
      + (a4*rho**4)/4 \
      + (a5*rho**5)/5 \
      + (a6*rho**6)/6 \
      + (a7*rho**7)/7 \
      + (a8*rho**8)/8 \
      + b1*G1 + b2*G2 + b3*G3 + b4*G4 + b5*G5 + b6*G6
    

    #eq.7 P_ = P*sigma**3/epsilon
    P_=rho*T \
     + a1*rho**2 + a2*rho**3 + a3*rho**4 + a4*rho**5 + a5*rho**6 + a6*rho**7 + a7*rho**8 \
     + F*(b1*rho**3 + b2*rho**5 + b3*rho**7 + b4*rho**9 + b5*rho**11)
    
    
    #Ur_ = Ur/N/epsilon
    #eq.9 
    Ur_=(c1*rho**1)/1 \
      + (c2*rho**2)/2 \
      + (c3*rho**3)/3 \
      + (c4*rho**4)/4 \
      + (c5*rho**5)/5 \
      + (c6*rho**6)/6 \
      + (c7*rho**7)/7 \
      + (c8*rho**8)/8 \
      + d1*G1 + d2*G2 + d3*G3 + d4*G4 + d5*G5 + d6*G6

    
    #eq.10
    Gr_ = Ar_ + P_/rho - T

    return Ar_, P_, Ur_, Gr_

class reduced_unit:
    #Usage: reduced_unit(300, 360, 2.71**3).cal_T() 
    def __init__(self,T,N,V,epsilon=93.20*8.314/1000, sigma=0.31589):
                #epsilon/k (K)   sigma (nm)
#SPC/E           78.20          0.31656
#tip4p           78.00          0.31540
#tip4p/2005      93.20          0.31589
#tip4p/ice      106.10          0.31668
        self.T=T
        self.N=N
        self.V=V
        self.epsilon=epsilon
        self.sigma=sigma
        
    def h(self):
        print("T-K, V-nm^3, epsilon-kJ/mol, sigma-nm")
    
    def cal_T(self):
        return self.T/(self.epsilon*1000/8.314)
    
    def cal_rho(self):
        return (self.N*self.sigma**3)/self.V
    
    def def_epsilon(self):
        return (self.epsilon)

    def def_sigma(self):
        return (self.sigma)

temp=float(sys.argv[1]) #K
Num=float(sys.argv[2]) 
boxsz=float(sys.argv[3])     #nm
epsilon=float(sys.argv[4])*8.314/1000
sigma=float(sys.argv[5])/10
deltaA=float(sys.argv[6])
#--------------------------------------------------
reduced_params=reduced_unit(temp, Num, boxsz**3, epsilon, sigma)
T_r=reduced_params.cal_T()
rho_r=reduced_params.cal_rho()
#---------------------------------------------------
result=MBWR(T_r, rho_r)

#--debug---
#print('P_reduced, U_res', result[1], result[2])
#----------

#---ideal term----
A_id =  np.log(360/(boxsz*10)**3) -1  #ln(rho*lambda**3) -1 
#----------------
print('A_res:' , result[0]/T_r )
print('A_id:'  , A_id )
print('A_LJ:'  , result[0]/T_r+A_id )
print('deltaA:', deltaA )
print('total:',deltaA+result[0]/T_r+A_id )
    
#T-K, N, V-nm^3, epsilon-kJ/mol, sigma-nm
EOF






    for i in `seq 0 1 $length`
    do
        if [ "${density}" != "no" ]; then 
            box_size=`echo "${N} ${density}"|awk '{printf "%.3f" ,($1*18/602/$2)^(1/3)}'`    
            Pressure[${i}]='none'
            ensemble=nvt
        fi
        echo "#sqrt(lambda) Coulomb energy" >>landa_coul_${Temperature[${i}]}_${Pressure[${i}]}.dat

        for f in `sed -n "2,$"p src/landa_weight.dat|awk '{print $1}'` #loop for sqrt(lambda)
        do
            coul=`echo "3 4 0"|gmx energy -f ${Temperature[${i}]}_${Pressure[${i}]}/${f}/rerun/rerun.edr -b 2000 -o ${Temperature[${i}]}_${Pressure[${i}]}/${f}/rerun/rerun.xvg|grep Coul|awk '{total+=$3} END {print total*1000/'${N}'/8.314/'${Temperature[${i}]}'}'`
            echo "${f} ${coul}" >>landa_coul_${Temperature[${i}]}_${Pressure[${i}]}.dat

        done
        deltaA=`paste landa_coul_${Temperature[${i}]}_${Pressure[${i}]}.dat src/landa_weight.dat|sed -n "2,$"p |awk '{print $2*$5}'|awk '{total+=$1} END {print total}'` 
        echo "#Free energy(NkT)---T=${Temperature[${i}]} K, P=${Pressure[${i}]} bar---" >>result.dat
        
        python3 MBWR.py ${Temperature[${i}]} ${N} ${box_size} ${epsilon} ${sigma} $deltaA >>result.dat
        echo "#------------------------------------------------">>result.dat
    done

fi

#Thermodynamic Integration
#"${Analyze} ="yes" &&  "${TI}" = "yes" 
#if [ "${TI}" = "yes" ]; then



#fi

