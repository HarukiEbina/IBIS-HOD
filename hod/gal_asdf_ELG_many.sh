#!/bin/bash -l
#SBATCH -J asdf_ELG
#SBATCH -N 1
#SBATCH -t 2:00:00
#SBATCH -o slurm_output/asdf.out
#SBATCH -e slurm_output/asdf.err
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH -A m68


module load python
conda activate test

export OMP_NUM_THREADS=8

# python gal_asdf_ELG_many.py 2.5
python gal_asdf_ELG_many.py 3.0