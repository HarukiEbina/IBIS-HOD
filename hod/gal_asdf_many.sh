#!/bin/bash -l
#SBATCH -J asdf
#SBATCH -N 1
#SBATCH -t 30:00
#SBATCH -o slurm_output/asdf.out
#SBATCH -e slurm_output/asdf.err
#SBATCH -q debug
#SBATCH -C cpu
#SBATCH -A m68


module load python
conda activate test

export OMP_NUM_THREADS=8

python gal_asdf_many.py 2.5
# python gal_asdf_many.py 3.0