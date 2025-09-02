#!/bin/bash -l
#SBATCH -J clustering
#SBATCH -N 1
#SBATCH -t 30:00
#SBATCH -o slurm_output/mm.out
#SBATCH -e slurm_output/mm.err
#SBATCH -q debug
#SBATCH -C cpu
#SBATCH -A desi


module load python
conda activate test

export OMP_NUM_THREADS=8

# python save_matter.py
# python ibis_clustering_mm.py 2.5
python ibis_clustering_mm.py 3.0

