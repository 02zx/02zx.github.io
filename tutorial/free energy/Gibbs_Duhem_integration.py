import subprocess
import time
import re
import os
import numpy as np


def generate_input_files(phase,temperature,pressure,seed):
    output_dir=f"{temperature}_{pressure}"
    os.makedirs(output_dir, exist_ok=True)
    path=f"{temperature}_{pressure}"
    in_file = f"{path}/in.{phase}"
    sh_file = f"{path}/{phase}.sh"
    with open(in_file, "w") as f:
        f.write("units           metal \n")
        f.write("dimension       3\n")
        f.write("boundary        p p p\n")
        f.write("atom_style      atomic\n")
        f.write(f"read_data       ../{phase}.data\n")
        f.write("labelmap atom 1 H 2 N 3 F\n")
        f.write("mass            1 1.008000\n")
        f.write("mass            2 14.00700\n")
        f.write("mass            3 18.99800\n")
        f.write("pair_style deepmd /home/Share/XZ/MLP/model-compress-pbe.pb\n")
        f.write("pair_coeff * * H N F\n")
        f.write("neighbor        2.0 bin\n")
        f.write("neigh_modify every 10 delay 0 check yes\n")
        f.write("thermo 10\n")
        f.write(f"dump rundump all custom 2000 ./{phase}.lammpstrj id type x y z\n")
        f.write("dump_modify rundump sort id\n")
        f.write(f"velocity all create {temperature} {seed}  dist gaussian\n")
        f.write(f'fix print_volume all print 10 "$(step) $(enthalpy) $(vol)" file thermo_{phase}.dat screen no title "# step Enthalpy(eV) Volume(A^3)" screen no\n')
        f.write(f"fix 1 all npt temp {temperature} {temperature} 0.1 aniso {pressure} {pressure} 1.0\n")
        f.write("timestep 0.0005\n")
        f.write("run 200000\n")
        f.write(f"write_data ./{phase}.data types labels\n")

    with open(sh_file, "w") as f:
        f.write('#!/bin/bash\n')
        f.write('#SBATCH -n 16\n')
        f.write('#SBATCH -N 1\n')
        f.write(f"#SBATCH -J {phase}_{temperature}_{pressure}\n")
        f.write('#SBATCH --gres=gpu:1\n')
        f.write('#SBATCH -p NVIDIAGeForceRTX3090\n')
        f.write('#SBATCH -e error.err\n')
        f.write('#SBATCH -o output.out\n')
        f.write('module purge\n')
        f.write('module add apps/lammps/2Aug2023-gpu3080\n')
        f.write(f"lmp -in ./in.{phase} -log {phase}.log")
        
    return path

def submit_jobs(path,phase):
    job_ids = []
    os.chdir(path)
    for inp in phase:
        cmd = ["sbatch", f"{inp}.sh"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        print(result)
        job_id = result.stdout.strip().split()[-1]  # "Submitted batch job 12345"
        job_ids.append(job_id)
    print(f"Submitted batch job: {job_ids}")
    os.chdir("../")
    return job_ids

def wait_for_jobs(job_ids):
    """轮询等待所有作业完成（用 squeue 检查）"""
    while True:
        cmd = ["squeue", "--jobs", ",".join(job_ids), "--noheader"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if not result.stdout.strip():  # 无输出表示所有作业完成
            print("###finished###")
            break
        time.sleep(30)  # 每 30 秒检查一次

def get_atom_num(path,filename):
    file_path=f"{path}/{filename}"
    with open(file_path, 'r') as f:
        content = f.read()
    
    # 正则匹配：开头可选空格 + 数字 + 空格 + atoms
    match = re.search(r'^\s*(\d+)\s+atoms', content, re.MULTILINE)
    atom_count = int(match.group(1))
    return atom_count
    
def slope_cal(path,phase):
    enthalpy=[]
    volume=[]
    for i in phase:
        data=np.loadtxt(f"{path}/thermo_{i}.dat")[10:]
        nmol=get_atom_num("./",f"{i}.data")/6
        H=np.mean(data[:,1])/nmol
        V=np.mean(data[:,2])/nmol
        enthalpy.append(H)
        volume.append(V)
        
    print(f"echo Properties phase 1 enthalpy {enthalpy[0]} volume {volume[0]} - phase 2 enthalpy {enthalpy[1]} volume {volume[1]}")
    conversion_factor_bar_eV_Acube=6.241509e-7
    slope=temperature*(volume[0]-volume[1])*conversion_factor_bar_eV_Acube/(enthalpy[0]-enthalpy[1])
    return slope
    

def flow(phase,temperature,pressure,seed):
    for p in phase:
        path=generate_input_files(p,temperature,pressure,seed)

    job_ids=submit_jobs(path,phase)
    wait_for_jobs(job_ids)
    print(f"echo kx calculation - temperature {temperature} - pressure {pressure}")
    slope=slope_cal(path,phase)
    return slope

def Runge_Kutta(phase,temperature,pressure,P_step,seed):
    k1=flow(phase,temperature,pressure,seed)
    k2=flow(phase,np.round(temperature+0.5*P_step*k1,2),pressure+0.5*P_step,seed)
    k3=flow(phase,np.round(temperature+0.5*P_step*k2),pressure+0.5*P_step,seed)
    k4=flow(phase,np.round(temperature+P_step*k3),pressure+P_step,seed)
    new_pressure=pressure+P_step
    new_temperature=np.round(temperature+(1./6.)*P_step*(k1+2*k2+2*k3+k4))
    return new_pressure,new_temperature


seed=79887
phase=["Ih","L"]
P_step=-1000

for loop in range(4):
    with open('thermal_condition.txt', 'r', encoding='utf-8') as f:
        lines = f.readlines()     
        if lines:                 
            last_line = lines[-1].strip()  

    pressure, temperature = map(float, last_line.split())

    P,T=Runge_Kutta(phase,temperature,pressure,P_step,seed)
    with open("thermal_condition.txt", "a") as f:
        f.write(f"{P} {T} \n")

