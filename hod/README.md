# Directory to run operations involving the entire AbacusHOD box

## Procedure to generate the necessary data products for analysis
### Note that if one only needs the final products, the existing .json files is sufficient. The only required step is to unzip the .json.gz files. 
1. Save the catalog of a large grid of HODs. This is necessary for Monte-Carlo fits in the `mc` directory.
   - `gal_asdf_many.sh`, `gal_asdf_ELG_many.sh` 
2. Compute the real-space auto-clustering of a large grid of HODs
   - `run_clustering_large_grid.sh`, `run_clustering_large_grid_ELG.sh`
3. Save the matter catalog of N-body boxes (`z = 2.5` and `3.0`). This is necessary for clustering measurements involving the matter field and to subtract the Poisson shot noise in PT fits.
   - `run_save_matter.sh`
4. Compute the matter-spectra for N-body
   - `run_clustering_mm.sh`
5. Compute auto- and cross-spectra (in real- and redshift-space, and with and without ZCV) for best-fit HOD models determined from Monte-Carlo fits. 
   - `run_clustering.sh`, `run_clustering_cross.sh`, `run_clustering_ELG.sh`, `run_clustering_cross_ELG.sh`

Code for perturbation theory fits to HOD clustering is in the `cobaya` subdirectory.
Key figures using the entire simulation box are generated in `hod_zcv.ipynb`.