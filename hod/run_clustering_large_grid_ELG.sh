#!/bin/bash -l
#SBATCH -J ELG2.5
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH -o clustering.out
#SBATCH -e clustering.err
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH -A desi


module load python
conda activate test

export OMP_NUM_THREADS=8


python ibis_clustering_large_grid_ELG.py '2.5'
# python ibis_clustering_large_grid_ELG.py '3.0'
