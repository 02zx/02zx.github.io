
```python
##input.py
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import numpy as np


pdb = PDBFile('/public/home/data_XZ/openMM/water.pdb')

forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HAngles)

platform = Platform.getPlatformByName('CUDA')
properties={'DeviceIndex': '0', 'Precision': 'single', 'UseBlockingSync':'false' }


integrator = NoseHooverIntegrator(300*kelvin, 5/picosecond, 0.002*picoseconds)
#integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

system.addForce(MonteCarloBarostat(pres*bar, 300*kelvin))
simulation = Simulation(pdb.topology, system, integrator, platform, properties)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
simulation.reporters.append(DCDReporter('output.dcd', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True, volume=True, speed=True))

simulation.step(1000000)
simulation.saveCheckpoint('state.chk')
```

```sh
#!/bin/bash
#SBATCH -n 1
#SBATCH -J openmm
#SBATCH -N 1
#SBATCH -p NVIDIAGeForceRTX3090
#SBATCH --gres=gpu:1
#SBATCH -e error.err1
#SBATCH -o output.out1


#/public/software/apps/openmm/8.1.1/lib/python3.12/site-packages/openmm/openmm.py
module add apps/openmm/8.1.1
#python3 -m openmm.testInstallation
python3 input.py

```
