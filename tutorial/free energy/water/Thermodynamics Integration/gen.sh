for i in `seq 225 10 443`
do
sed -i s"@ref_t=.*@ref_t=$i@" eq.mdp

gmx_mpi grompp -f eq.mdp -c conf.gro -p spce.top -o $i/md

sed -i s"@i=.*@i=$i@" run.sh
sbatch run.sh
done



