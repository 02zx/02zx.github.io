NepTrain
```bash
NepTrain perturb train.xyz -n 20000 -c 0.03 -d 0.2 #generating pertubation data
NepTrain select perturb.xyz -max 100 -d 0.1 #far point sampling
#The number of sampled confs determined by -d
```
run NepTrain using the settings in job.yaml
```bash
NepTrain train job.yaml
```
the directory contains:
cache/  dpdispatcher.log  INCAR  job.yaml  restart.yaml  structure/  train.xyz

structure directory contains the initial configuration for MD sampling.

cache contains the output of active learning

Initial train.xyz is the training set

INCAR is the input of VASP calculation.

job.yaml
```yaml
version: 2.2.0
dft_job: 10 #The number of tasks submitted when calculating single-point energy with DFT.
#所有任务提交的根目录
gpumd_split_job: temperature  #Split gpumd tasks by temperature or structure
work_path: ./cache  #Root directory for all task submissions.
current_job: nep
#If the current_job has three states: nep, gpumd, dft, and if train.xyz has not been calculated, set it to vasp; otherwise, directly set it to use nep to train the potential function, or use gpumd.
#从nep训练开始->gpumd->select->dft->nep...
#
generation: 1  #Marking resume tasks.
init_train_xyz: ./train.xyz  #Initial training set; if not calculated, set current_job to vasp.
init_nep_txt: ./nep.txt  #If current_job is set to gpumd, a potential function must be provided; otherwise, it can be ignored.
nep:
  #Does it support restarting? If true, the potential function for the next step will continue from this step for nep_restart_step steps.
  #The program will automatically set lambda_1 to 0.
  #If false, retrain from scratch every time.
  nep_restart: true
  nep_restart_step: 20000
  #Optional; if you need to modify the number of steps, simply provide a file in the current path.
  #If there is no such file, the number of steps will be automatically generated based on the training set.
  nep_in_path: ./nep.in
  #Optional
  test_xyz_path: test.xyz
  machine:
    #https://docs.deepmodeling.com/projects/dpdispatcher/en/latest/context.html
    context_type: LazyLocal
    #https://docs.deepmodeling.com/projects/dpdispatcher/en/latest/batch.html
    batch_type: Slurm
    local_root: ./
    remote_root: ~/neptrain/
    remote_profile:
      hostname: ''
      username: ''
      key_filename: ''
      port: 22


  resources:
    number_node: 1
    cpu_per_node: 32
    gpu_per_node: 2
    queue_name: NVIDIAGeForceRTX3090 
    group_size: 1
    custom_flags:
    - '#SBATCH --job-name=NepTrain-NEP'

    prepend_script:
    - ' '



dft:
  software: vasp #you can switch vasp or abacus
  cpu_core: 32
  kpoints_use_gamma: true  #ASE defaults to using M-point k-mesh, but here we default to using the gamma-centered grid; this can be set to false.

  incar_path: auto  #A path should be passed in. If it is auto, the corresponding file name will be switched according to the software, such as INCAR or INPUT

  use_k_stype: kspacing
  #--ka
  kpoints:
  - 20   #a
  - 20   #b
  - 20   #c
  kspacing: 0.2
  machine:
    context_type: LazyLocal
    batch_type: Slurm
    local_root: ./
    remote_root: ~/neptrain/
    remote_profile:
      hostname: ''
      username: ''
      key_filename: ''
      port: 22


  resources:
    number_node: 1
    cpu_per_node: 64
    gpu_per_node: 0
    queue_name: INTEL_8358P
    group_size: 1
    custom_flags:
    - '#SBATCH --job-name=NepTrain-dft'

    prepend_script:
    - ' '


gpumd:
#Time for iterative progressive learning in units of picoseconds.
#The first active learning is at 10ps, the second at 100ps, with a total of four active learning sessions.
  step_times:
  - 10
  - 100
  - 500
  - 1000
#Each time active learning is performed, all structures in model_path will undergo molecular dynamics (MD) simulations at the following temperatures, followed by sampling.
  temperature_every_step:
  - 50
  - 100
  - 150
  - 200
  - 250
  - 300
  model_path: ./structure
  run_in_path: ./run.in

  machine:
    context_type: LazyLocal
    batch_type: Slurm
    local_root: ./
    remote_root: ~/neptrain/
    remote_profile:
      hostname: ''
      username: ''
      key_filename: ''
      port: 22


  resources:
    number_node: 1
    gpu_per_node: 1
    queue_name: NVIDIAGeForceRTX3090
    group_size: 1
    custom_flags:
    - '#SBATCH --cpus-per-task=16'
    - '#SBATCH --job-name=NepTrain-gpumd'
    prepend_script:
    - ' '

select:
  #After completing this round of MD, a maximum of max_selected structures will be selected from all trajectories.
  max_selected: 50
  min_distance: 0.01   #Hyperparameters for farthest point sampling
  filter: 0.6    #Passing a coefficient enables bond length detection, and bonds shorter than the sum of the covalent radii multiplied by the coefficient are considered unphysical structures.
  machine:
    context_type: LazyLocal
    batch_type: Slurm
    local_root: ./
    remote_root: ~/neptrain/
    remote_profile:
      hostname: ''
      username: ''
      key_filename: ''
      port: 22


  resources:
    #Select mainly uses CPU
    number_node: 1
    gpu_per_node: 0
    queue_name: INTEL_8358P
    group_size: 1
    custom_flags:
    - '#SBATCH --cpus-per-task=1'
    - '#SBATCH --job-name=NepTrain-select'

    prepend_script:
    - ' '


limit:
  force: 20  #Limit the force of the structure to between -force and force
```

GPUMD inputs
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
