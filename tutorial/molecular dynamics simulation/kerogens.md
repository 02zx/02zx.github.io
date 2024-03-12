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



## 计算helium void fraction

## 计算ideal gas Rosenbluth weight

p=0-24MPa, T=338K

[^1]:页岩气吸附与CO_2驱替及封存机理的分子模拟研究-周娟
[^2]:https://www.materialsdesign.com/Publications/Ungerer2015
[^3]:Wang, C.; Li, W.; Liao, K.; Wang, Z.; Wang, Y.; Gong, K. AuToFF Program, Vesrion 1.0. Hzwtech. Shanghai 2023, see https://cloud.hzwtech.com/web/product-service?id=36.
