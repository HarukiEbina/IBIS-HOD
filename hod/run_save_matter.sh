#!/bin/bash -l
#SBATCH -J clustering
#SBATCH -N 1
#SBATCH -t 4:00:00
#SBATCH -o slurm_output/save_matter.out
#SBATCH -e slurm_output/save_matter.err
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH -A desi


module load python
conda activate test

export OMP_NUM_THREADS=8

# python save_matter.py 2.5
python save_matter.py 3.0
# python ibis_clustering_cross.py
