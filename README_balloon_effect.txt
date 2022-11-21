21 November 2022

Direct all requests to Eleonora Dávalos edavalosa@eafit.edu.co

These scripts and data reproduce the results of Diffusion of crime control benefits: Forced eradication and coca crops in Colombia. 

Required STATA packages:

shp2dta
estout
xsmle



1) Place all scripts in the same folder as the data.

2) Decompress files

3) Stata do file 1. GenMatrices_0522: This do file produce and export all matrices used in spatial regressions

4) Stata do file 2. Estimations_0522: This do files produce summary regression results and wald statistics Table 4, Table 5, Table 6, table 7

5) Matlab m file 3. MatlabEstimations: This m file compute bayesian probabilities presented in Table 5 and Table 6. 
   This code is based on Elhorst JP (2010) Spatial Panel Data Models, and the scripts provided https://spatial-panels.com/ and Paul Elhorst website.
   All auxiliary functions can be found in https://spatial-panels.com/.

6) Matlab m file 4. LMTestPaperIlicitCrops2022: This m file computes lagrange multiplier test presented in Table 4.
   This code is based on Elhorst JP (2010) Spatial Panel Data Models, and the scripts provided https://spatial-panels.com/ and Paul Elhorst website.
   All auxiliary functions can be found in https://spatial-panels.com/.
