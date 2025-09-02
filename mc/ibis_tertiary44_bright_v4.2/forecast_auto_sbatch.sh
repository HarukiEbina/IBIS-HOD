#!/bin/bash -l
#SBATCH -J forecast_sbatch
#SBATCH -N 1
#SBATCH -t 00:15:00
#SBATCH -o clustering.out
#SBATCH -e clustering.err
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH -A m68

module load python
conda activate myenv_abacus

cd /pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/mc

# argv[1] is path to config file
date
echo /pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/mc/AbacusSummit_high_c000_ph100/z3.0/ibis_tertiary44_bright_v4.2/logM_cut_11.00_logM1_12.00_sigma_0.66_kappa_1.00_alpha_0.66/configs.json
python save_hod.py /pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/mc/AbacusSummit_high_c000_ph100/z3.0/ibis_tertiary44_bright_v4.2/logM_cut_11.00_logM1_12.00_sigma_0.66_kappa_1.00_alpha_0.66/configs.json
date
