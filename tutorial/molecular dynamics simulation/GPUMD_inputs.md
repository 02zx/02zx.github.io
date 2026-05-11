```bash
#!/bin/bash
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -J Ih_nep_scan_0.0
#SBATCH --gres=gpu:1
#SBATCH -p NVIDIAGeForceRTX3090
#SBATCH -e error.err
#SBATCH -o output.out


module load compiler/gcc/7.3.1
module load compiler/intel/2021.3.0
module load mpi/intelmpi/2021.3.0
module load mathlib/fftw/intelmpi/3.3.9_single

conf=conf.xyz
T=500
p=0.0
run=200000 #0.1ns

mkdir ${T}_${p}
cat > ./${T}_${p}/run.in<<EOF
potential   nep.txt
velocity    ${T}
time_step   0.5
ensemble    npt_scr ${T} ${T} 100 ${p} ${p} ${p} 100 100 100 2000
#the stochastic cell rescaling T_init T_end tau_T(step) p_xx(GPa) p_yy p_zz C_xx C_yy C_zz tau_p(step)
dump_exyz   20000
#dump_exyz <interval> <has_velocity> <has_force> <has_potential> <separated>
dump_thermo 2000
#column   1 2 3 4  5  6  7   8   9   10 11 12 13 14 15 16 17 18
#quantity T K U Pxx Pyy Pzz Pyz Pxz Pxy ax ay az bx by bz cx cy cz
run         ${run}
EOF

cp ${conf} ./${T}_${p}/model.xyz
cd ${T}_${p}
/public/home/data_XZ/GPUMD/GPUMD-4.7/src/gpumd


```
