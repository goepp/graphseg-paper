#!/bin/sh
#SBATCH -J lauchRscript
#SBATCH --job-name=Infer5
#SBATCH --mail-user=vivien.goepp@mines-paristech.fr
#SBATCH --mail-type=END,FAIL
#SBATCH -o infer5.out
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --time=5-23:59:59
#SBATCH --mem=30gb
pwd; hostname; date
 
echo "Running on $SLURM_CPUS_ON_NODE CPU cores"
#Rscript infer_any.R neth municip province
#Rscript infer_any.R utrecht district municip
#Rscript infer_any.R utrecht neigh province
#Rscript infer_any.R utrecht neigh municip
Rscript infer_any.R utrecht neigh district
#Rscript infer_any.R neth neigh municip

