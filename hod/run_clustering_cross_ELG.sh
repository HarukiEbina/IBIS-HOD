#!/bin/bash -l
#SBATCH -J gm
#SBATCH -N 1
#SBATCH -t 30:00
#SBATCH -o slurm_output/clustering.out
#SBATCH -e slurm_output/clustering.err
#SBATCH -q debug
#SBATCH -C cpu
#SBATCH -A desi


module load python
conda activate test

export OMP_NUM_THREADS=8

# python save_matter.py
python ibis_clustering_cross_ELG.py
