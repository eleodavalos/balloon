21 November 2022

Direct all requests to Eleonora Dávalos edavalosa@eafit.edu.co

These scripts and data reproduce the results of Diffusion of crime control benefits: Forced eradication and coca crops in Colombia. 


EDITAR DE AQUÍ PARA ABBAJO

Required STATA packages:

reghdfe
estout
movestay
svmat2


1) Place all scripts in the same folder as the data.

2) sem_dat_v7.R tabulates the data and includes the adjacency matrix for downstream analyses, saves a RData file.

3) piecewise_explore_v8.R generates piecewise SEM,

4) sem_bet_v8nb.R and sem_bet_v8zi.R run Bayesian models, these run parallel chains and take <7 days in a cluster but 4 months on a mini-cluster.

5) conditional_v4 and smooth_v4.R process the Bayesian results and save their corresponding files.

6) plot_summaries_v9.R, plot_smooth_v4.R, plot_conditional_v5.R use already generated RData files to plot coefficients, GMRF smooths, and conditional effects.

7) loo_pit_plo_v9.R generates and plots loo-pit simulations and requires both models to be compared.

8) network1_v3.R makes a network plot based on cov_brms_v3.csv, which was compiled by hand. The colors were adjusted by hand.





1) Place all scripts in the same folder as the data.

2) Decompress files

3) Descriptive tables_Final_JAAE: this do file produces the descriptive statistics of agricultural units presented in Table 1 and Table 2 of the paper.

4) First Stage, OLS and IV regressions_Final_JAAE: this do file produce the First Stage and Second Stage regressions presented in Table 3 and Table 4 of the paper.

5) MTE regression_Final_JAAE: this do file produce the MTE estimation presented in Table 5 and Figure 1 of the paper.

6) Heterogenous effects by land size_Final_JAAE: this do file produce results of MTE by land size presented in Table 6.



