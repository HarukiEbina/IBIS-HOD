#!/bin/bash -l
#SBATCH -J gg
#SBATCH -N 1
#SBATCH -t 30:00
#SBATCH -o slurm_output/zcv.out
#SBATCH -e slurm_output/zcv.err
#SBATCH -q debug
#SBATCH -C cpu
#SBATCH -A desi


module load python
conda activate test

export OMP_NUM_THREADS=8

python ibis_clustering_ELG.py