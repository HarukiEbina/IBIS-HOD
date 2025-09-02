#!/bin/bash -l
#SBATCH -J sim_min
#SBATCH -t 10:00
#SBATCH -N 1
#SBATCH -q debug
#SBATCH -C cpu
#SBATCH -A desi
#SBATCH -o slurm_output/fs_min_sim.txt

date
#
module load python
source activate cobaya_07292024

# export COBAYA_USE_FILE_LOCKING=False
#
# export PYTHONPATH=/global/homes/h/hebina/cobaya_runner:$PYTHONPATH
# export PYTHONPATH=/global/homes/h/hebina/cobaya_runner/likelihoods:$PYTHONPATH
export PYTHONPATH=/pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/hod/cobaya:$PYTHONPATH
export PYTHONPATH=/pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/hod/cobaya/likelihoods:$PYTHONPATH
# export PYTHONPATH=${PYTHONPATH}:/pscratch/sd/m/mmaus/ShapeFit_DESI_Cutsky/emulator/fullshape/abacus_fid/wcdm/redef
# # export PYTHONPATH=${PYTHONPATH}:/global/homes/m/mmaus/Python/velocileptors
# export PYTHONPATH=${PYTHONPATH}:/pscratch/sd/m/mmaus/Cobaya/Packages/code
export OMP_NUM_THREADS=8
#
# rm ./output/official_w0wa/Ib_Y7_scen2_PR3.input.yaml.locked

echo "Setup done.  Starting to run code ..."

srun -n 16 -c $OMP_NUM_THREADS cobaya-run -f yamls/fs_minimize_z2.5.yaml

# #
date
#
