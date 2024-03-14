# II-C型干酪根体系建模[^1]

分子取自[^2], 力场:GAFF2+AM1-BCC获取自autoFF[^3]

输入文件在kerogens.zip中.

## 模拟流程:

1.对单分子进行能量最小化(最速下降法)
```bash
mkdir em eq prod

#---修改topo.top文件----
[ molecules ]
; Molecule      nmols
IIC        1
#------------

gmx grompp -f em.mdp -c IIC.pdb -p topo.top -o em/em
gmx mdrun -v -deffnm em/em
```

2.在300K, 1bar下进行npt平衡400ps
```bash
gmx grompp -f eq.mdp -c em/em.gro -p topo.top -o eq/eq
sbatch run.sh 
```

3.在10 $^3$ nm盒子内插入8个分子
```bash
gmx insert-molecules -box 10 10 10 -ci eq/eq.gro -nmol 8 -o conf.pdb

#---修改topo.top文件----
[ molecules ]
; Molecule      nmols
IIC        8
#------------

gmx grompp -f eq.mdp -c conf.pdb -p topo.top -o prod/eq
sbatch run.sh
```

4.在300K, 200bar(20MPa)下进行npt平衡400ps
```bash
#对eq.mdp进行修改
#----eq.mdp------------------------
Pcoupl              =  berendsen
Pcoupltype          =  isotropic
tau_p               =  1
compressibility     =  4.5e-5
ref_p= 200 ;bar 
#--------------------------------

gmx grompp -f eq.mdp -c prod/eq.gro -p topo.top -o prod/eq2
sbatch run.sh
```

5.将体系升温到900K后退火到338.15K
```bash
gmx grompp -f prod.mdp -c prod/eq2.gro -p topo.top -o prod/md
sbatch run.sh
```

# 含水量




# 甲烷吸附

需要准备的文件:
```bash
kerogens.cif
pseudo_atoms.def
force_field_mixing_rules.def
helium.def
ch4.def
```
---kerogens.cif---

将上一步中模拟产生的prod/md.gro转换为pdb文件, 然后再用Material Studio打开保存成cif即可.
文件见附件

---pseudo_atoms.def---

GAFF2_XX的参数取自AutoFF, 下同

```bash
#number of pseudo atoms. Created by AutoFF
2
# type print as scat oxidation mass charge polarization B-factor radii connectivity anisotropic anisotropic-type tinker-type
 Tra_CH4 yes    CH4 CH4 0  16.043000    0.000000 0.0 1.0 1.0 0 0.0 relative 0
      He yes    He  He  0   4.002602    0.0      0.0 1.0 1.0 0 0   relative 0
```
---force_field_mixing_rules.def---

```bash
# general rule for shifted vs truncated. Created by AutoFF
shift
# general rule tail corrections
no
# number of defined interactions
17
# type interaction, parameters. IMPORTANT: define generic matches first
GAFF2_ca  lennard-jones    49.718090    3.315212
GAFF2_c2  lennard-jones    49.718090    3.315212
GAFF2_na  lennard-jones   102.757429    3.205810
GAFF2_hn  lennard-jones     5.032195    1.106496
GAFF2_ha  lennard-jones     8.101834    2.625479
GAFF2_h4  lennard-jones     8.101834    2.536389
GAFF2_c3  lennard-jones    54.247066    3.397710
GAFF2_hc  lennard-jones    10.466966    2.600177
GAFF2_nb  lennard-jones    47.352958    3.384168
GAFF2_os  lennard-jones    36.533738    3.156098
 GAFF2_c  lennard-jones    49.718090    3.315212
GAFF2_ss  lennard-jones   142.109244    3.532413
 GAFF2_o  lennard-jones    73.621018    3.048121
GAFF2_oh  lennard-jones    46.799417    3.242871
GAFF2_ho  lennard-jones     2.365132    0.537925
 Tra_CH4  lennard-jones   147.999944    3.730000
      He  lennard-jones    10.9         2.64
# general mixing rule for Lennard-Jones
Lorentz-Berthelot
```


---helium.def---

helium.def和ch4.def可以直接从RASPA的安装路径下获取.

```bash
# critical constants: Temperature [T], Pressure [Pa], and Acentric factor [-]
5.2
228000.0
-0.39
# Number Of Atoms
1
# Number Of Groups
1
# Alkane-group
flexible
# number of atoms
1
# atomic positions
0 He
# Chiral centers Bond  BondDipoles Bend  UrayBradley InvBend  Torsion Imp. Torsion Bond/Bond Stretch/Bend Bend/Bend Stretch/Torsion Bend/Torsion IntraVDW IntraCoulomb
               0    0            0    0            0       0        0            0         0            0         0               0            0        0            0
# Number of config moves
0
```
---ch4.def---
```bash
# critical constants: Temperature [T], Pressure [Pa], and Acentric factor [-]. Created by AutoFF
0.0
0.0
0.0
#Number Of Atoms
1
#Number Of Group
1
# ch4-group
flexible
#number of atoms
1
# atomic positions
0      Tra_CH4      1.064000   -0.028700    0.071300
# Chiralcenters Bond BondDipoles Bend UrayBradley InvBend Torsion Imp.Torsion Bond/Bond Stretch/Bend Bend/Bend Stretch/Torsion Bend/Torsion IntraVDW IntraCoulomb
0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Number of config moves
0
```

## 计算干酪根体系的helium void fraction
输入文件内容如下(参考自手册example):

由于unitcell体积比较大( $3^3$ nm $^3$ ), 模拟会比较耗时

---He_frac.in---
```bash
SimulationType		MonteCarlo
NumberOfCycles		500000
printEvery		10000
PrintPropertiesEvery	10000

Forcefield		local

Framework		0
FrameworkName		kerogens
UnitCells		1 1 1
ExternalTemperature	338.0

Component 0 MoleculeName		helium
	    MoleculeDefinition		local
	    WidomProbability		1.0
	    CreateNumberOfMolecules 	0
```



## 计算甲烷的ideal gas Rosenbluth weight

p=0-24MPa, T=338K

[^1]:页岩气吸附与CO_2驱替及封存机理的分子模拟研究-周娟
[^2]:https://www.materialsdesign.com/Publications/Ungerer2015
[^3]:Wang, C.; Li, W.; Liao, K.; Wang, Z.; Wang, Y.; Gong, K. AuToFF Program, Vesrion 1.0. Hzwtech. Shanghai 2023, see https://cloud.hzwtech.com/web/product-service?id=36.
