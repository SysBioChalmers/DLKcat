#!/bin/bash
#SBATCH -A C3SE2021-1-16
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -o out.txt
#SBATCH -t 10:00:00
#SBATCH --mail-user=feiranl@chalmers.se
#SBATCH --mail-type=end
module load GCCcore/8.3.0
module load MATLAB intel/2018b GMP
module load Gurobi/8.0.0
a1=1
b1=40
matlab -nodesktop -singleCompThread -r "writeFile_pre_cluster($a1,$b1)"
