# LBODE Extended Model
Code for Multi-Cellular Network Model Predicts Alterations in Glomerular Endothelial Structure in Diabetic Kidney Disease

[![DOI]()]()

## Overview
This extended logic-based ODE model predicts the effects of long-term exposure to glucose on pro-inflammatory macrophages and glomerular endothelial cells in diabetic mice. A novel model-based approach to predict the effect of glucose and disease intervention strategies on changes in fenestration structure in glomerular endothelial cells in diabetic kidney disease.

## Authors
Krutika Patidar<sup>a</sup>,  Ashlee N. Ford Versypt<sup>a,b,c,d</sup>

<sup>a</sup>Department of Chemical and Biological Engineering, University at Buffalo, The State University of New York, Buffalo, NY, USA<br/>
<sup>b</sup>Department of Biomedical Engineering, University at Buffalo, The State University of New York, Buffalo, NY, USA<br/>
<sup>c</sup>Institute for Artificial Intelligence and Data Science, University at Buffalo, The State University of New York, Buffalo, NY, USA<br/>
<sup>d</sup>Department of Pharmaceutical Sciences, University at Buffalo, The State University of New York, Buffalo, NY, 14215, USA<br/>

## Manuscript
K. Patidar and A. N. Ford Versypt, Multi-Cellular Network Model Predicts Alterations in Glomerular Endothelial Structure in Diabetic Kidney Disease, bioRxiv preprint, DOI: 10.1101/2024.12.30.630833 [Preprint](https://www.biorxiv.org/content/10.1101/2024.12.30.630833v1)

## Scripts

* call_ODE_IVV_model.m calls LBODE extended model and performs user-specific task.
* coupledODE_IVV_step.m provides a list of LBODEs for the extended model.
* coupledODE_physParams.m consists a set of parameters that define change in fenestration structure.
* coupledODE_IVIVC_params.m consists a set of extended LBODE model parameters.
* coupledODE_IVV_run.m simulates the extended LBODE model for mean glucose input without variability.
* step_function.m simulates a step change in uniformly sampled glucose concentration between 6-20 weeks.
* coupledODE_IVV_multirun.m simulates the extended LBODE model for variability in glucose input and variability in estimated parameters.
* multistart_param_opt.m is a multi-start non-linear least squares optimization routine for parameter estimation.
* coupledODE_IVV_SSE.m computes the mean sum of squared error (SSE) between predicted output and observed data.
* post_IVV.m simulates in silico predictions post optimization.
* pub_plots.m recreates the figures in the publication/pre-print.
* LHS_Call.m This supporting function uses Latin hypercube sampling to create sample subsets of parameters within a given range.
* jbfill.m This supporting function will fill a region with a color between the two vectors.
* sig_star.m is a script to add statistical significance stars to bar charts.


## Data
data folder consists of mat and csv files used to simulate or optimize the model.

* finch_density_ctrl.mat consists data for fenestration density for healthy mice.
* finch_fenestration_disease.mat consists data for fenestration density and width for diabetic mice.
* GLU_data.mat consists lumped data for observed glucose concentration in diabetic mice from multiple studies.
* dbmice_density_population.csv consists data for fenestration density in individual diabetic mice subject.
* dbmice_diameter_population.csv consists data for fenestration diameter in individual diabetic mice subject. 
* MC_25_fen_multirun.mat provides saved model variables and output to recreate 95% credible intervals from prediction posteriors.
 
## Recommended Supplementary Packages
* [Netflux](https://github.com/saucermanlab/Netflux) is a package that generates certain equations and utility functions for coupledODE_IVV_step.m

## Acknowledgements
Research reported in this publication was supported by the National Institute of General Medical Sciences of the National Institutes of Health under award number R35GM133763 and NSF CAREER
2133411. The content is solely the responsibility of the authors and does not necessarily represent the official views of the funding agencies.
