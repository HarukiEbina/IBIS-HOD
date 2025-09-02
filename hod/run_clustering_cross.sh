#!/bin/bash -l
#SBATCH -J gm
#SBATCH -N 1
#SBATCH -t 30:00
#SBATCH -o slurm_output/cross.out
#SBATCH -e slurm_output/cross.err
#SBATCH -q debug
#SBATCH -C cpu
#SBATCH -A m68


module load python
conda activate test

export OMP_NUM_THREADS=8

# python save_matter.py
python ibis_clustering_cross.py 3.0
python ibis_clustering_cross.py 2.5
