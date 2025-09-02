
module load python
conda activate myenv_abacus

cd /pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/mc
export OMP_NUM_THREADS=16
# argv[1] is path to config file
echo /pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/mc/AbacusSummit_high_c000_ph100/z3.0/ibis_tertiary44_clustering_v4.2/logM_cut_11.75_logM1_12.75_sigma_0.33_kappa_1.00_alpha_0.33_Q_100.00_gamma_1.00_p_max_0.66/configs.json
python save_hod.py /pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/mc/AbacusSummit_high_c000_ph100/z3.0/ibis_tertiary44_clustering_v4.2/logM_cut_11.75_logM1_12.75_sigma_0.33_kappa_1.00_alpha_0.33_Q_100.00_gamma_1.00_p_max_0.66/configs.json
            