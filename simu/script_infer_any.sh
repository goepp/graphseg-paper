#!/bin/sh
#SBATCH -J lauchRscript
#SBATCH --job-name=SpatialInferAny
#SBATCH --mail-user=vivien.goepp@mines-paristech.fr
#SBATCH --mail-type=END
#SBATCH -o infer_any.out
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --time=23:59:59
#SBATCH --mem=2gb
pwd; hostname; date
 
echo "Running on $SLURM_CPUS_ON_NODE CPU cores"
Rscript infer_any.R

