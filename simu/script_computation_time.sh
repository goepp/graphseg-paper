#!/bin/sh
#SBATCH -J lauchRscript
#SBATCH --job-name=SpaCompTime
#SBATCH --mail-user=vivien.goepp@mines-paristech.fr
#SBATCH --mail-type=END,FAIL
#SBATCH -o computation_time.out
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --time=10-23:59:59
#SBATCH --mem=4gb
pwd; hostname; date
 
echo "Running on $SLURM_CPUS_ON_NODE CPU cores"
Rscript computation_time.R

