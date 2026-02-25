
## Input for lammps
---input.in---
```bash
clear
# --------------------- VARIABLES-------------------------
variable        NSTEPS          equal 1000000000        #not important
variable        THERMO_FREQ     equal 10
variable        DUMP_FREQ       equal 10
variable        NREPEAT         equal ${NSTEPS}/${DUMP_FREQ}
#variable        TEMP            equal 100
#variable        PRES            equal 10*10000          #5GPa
#variable        TAU_T           equal 0.100000
#variable        TAU_P           equal 0.500000
#variable        DENS            equal 432*18/6.02/vol*10000 #kg/m3
# ---------------------- INITIALIZAITION ------------------
units           metal
boundary        p p p
atom_style      atomic
# --------------------- ATOM DEFINITION ------------------
#box             tilt large
read_data ./confs/III.data
labelmap atom 1 H 2 N 3 F
mass 1 1.008000 #H
mass 2 14.007  #N
mass 3 18.998   #F
pair_style deepmd /home/Share/XZ/MLP/model-compress-pbe.pb
pair_coeff * * H N F
neighbor        2.0 bin
neigh_modify every 10 delay 0 check yes

timestep        0.0002
fix             1 all ipi VARADDRESS 32345 unix
# --------------------- RUN ------------------------------
run             ${NSTEPS}
#write_data      npt.data
```

## Input for i-Pi
---input.xml---
```bash
<simulation verbosity='high'>

  <output prefix='sim'>
    <properties stride='10' filename='therm'>
     [step, time{picosecond}, temperature{kelvin},
      pressure_md{gigapascal}, pressure_cv{gigapascal},
      volume{angstrom3}, cell_abcABC] </properties>
    <properties stride='10' filename='ene'>
     [step, conserved{electronvolt},
       potential{electronvolt}, kinetic_cv{electronvolt}]  </properties>
    <trajectory filename='pos' format='pdb' stride='100' cell_units='angstrom'> positions{angstrom} </trajectory>
    <trajectory filename='pos-ct' format='pdb' stride='100' cell_units='angstrom'> x_centroid{angstrom} </trajectory>
    <trajectory filename='rg' format='pdb' stride='100' cell_units='angstrom'> r_gyration{angstrom} </trajectory>
    <checkpoint stride='200'/>
  </output>

  <total_steps>1000000000</total_steps>

  <prng>
    <seed>32345</seed>
  </prng>

  <ffsocket name='lammps' mode='unix'>
    <address> VARADDRESS </address>
  </ffsocket>

  <system>
    <initialize nbeads='32'>
      <file mode='pdb'> confs/eq.pdb </file>
      <velocities mode='thermal' units='kelvin'> 500 </velocities>
    </initialize>
    <forces>
      <force forcefield='lammps'> </force>
    </forces>
    <motion mode='dynamics'>
      <dynamics mode='npt'>
        <barostat mode="flexible">
          <tau units="femtosecond"> 200 </tau>
          <thermostat mode="langevin">
            <tau units="femtosecond"> 100 </tau>
          </thermostat>
        </barostat>
        <timestep units='femtosecond'> 0.2 </timestep>
        <thermostat mode='pile_g'>
          <tau units='femtosecond'> 10 </tau>
        </thermostat>
      </dynamics>
    </motion>
    <ensemble>
      <temperature units='kelvin'> 500 </temperature>
      <pressure units="gigapascal"> 2.0 </pressure>
    </ensemble>
  </system>

</simulation>
```

## script for simulation
---run.sh---
```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --gres=gpu:2
#SBATCH -p NVIDIAGeForceRTX3090
#SBATCH -J 2.0GPa_500K_PIMD_III
#SBATCH -e error.err
#SBATCH -o output.out

module purge
module load apps/lammps/2Aug2023-gpu-dp-ipi

# ======= ipi running =======

rm /tmp/ipi_VARPRES-VARTEMP
i-pi input.xml &> log.i-pi &

sleep 5

rm -fr log-lmp
mkdir log-lmp

echo "Launching 8 clients on GPU 0..."
for nlmp in `seq 0 7`; do
    export CUDA_VISIBLE_DEVICES=0
    lmp -in npt.in &> log-lmp/log.lmp.$nlmp &
done

echo "Launching 8 clients on GPU 1..."
for nlmp in `seq 8 15`; do
    export CUDA_VISIBLE_DEVICES=1
    lmp -in npt.in &> log-lmp/log.lmp.$nlmp &
done

wait

# ======= ipi job done =======
```
