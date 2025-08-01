#!/bin/bash -l
#SBATCH -J large_grid
#SBATCH -N 1
#SBATCH -t 6:00:00
#SBATCH -o slurm_output/large_grid.out
#SBATCH -e slurm_output/large_grid.err
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH -A desi


module load python
conda activate test

export OMP_NUM_THREADS=8

python ibis_clustering_large_grid.py 2.5
# python ibis_clustering_large_grid.py 3.0
# python ibis_clustering_large_grid_velbias.py