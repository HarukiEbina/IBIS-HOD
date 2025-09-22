# Directory for Monte-Carlo fits to data

Run `make_many_forecasts.py` in each directory to make many realizations of the data using the AbacusSummit box and compute the mean and covariance of clustering. 
This step is computationally expensive, so the first recommendation would be to use the existing files `r0_arr.json` and `r0_arr_ELG.json` for both the $r_0$ and $\chi^2$ information of a large grid of simulations. 
`bestfit.json` points to the best-fit HOD model identified in the analysis.

The goodness-of-fit of these to data are calculated in corresponding ipython notebooks, but this requires the simulations to be run, which are computationally expensive. 
One can instead arrive at the final product, e.g.~a $r_0$ vs.~$\chi^2$ distribution in Figure 7 and 15 of Ebina et al.~2025, by running the 

