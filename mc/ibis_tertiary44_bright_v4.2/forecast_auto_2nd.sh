
module load python
conda activate myenv_abacus

cd /pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/mc
export OMP_NUM_THREADS=16
# argv[1] is path to config file
echo /pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/mc/AbacusSummit_high_c000_ph100/z3.0/ibis_tertiary44_bright_v4.2/logM_cut_12.00_logM1_14.00_sigma_0.66_kappa_1.00_alpha_0.66_Q_100.00_gamma_5.00_p_max_0.66/configs.json
python save_hod.py /pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/mc/AbacusSummit_high_c000_ph100/z3.0/ibis_tertiary44_bright_v4.2/logM_cut_12.00_logM1_14.00_sigma_0.66_kappa_1.00_alpha_0.66_Q_100.00_gamma_5.00_p_max_0.66/configs.json
            