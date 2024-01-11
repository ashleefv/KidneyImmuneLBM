# KidneyImmuneLBM
Code for Logic-Based Modeling of Inflammatory Macrophage Crosstalk with Glomerular Endothelial Cells in Diabetic Kidney Disease

[![DOI](https://zenodo.org/badge/642046465.svg)](https://zenodo.org/badge/latestdoi/642046465)

## Overview
This logic-based ODE model predicts the effects of glucose and inflammatory stimulus on pro-inflammatory macrophages and glomerular endothelial cells in diabetic kidney disease. A protein signaling network describes the crosstalk between macrophages and glomerular endothelial cells stimulated by glucose and LPS, and it consists of 30 species and 40 interactions. The model inputs (glucose or LPS) are 0 or 1 when the input is inactive or fully active. The model species hold a value between 0 and 1. A set of 30 differential equations define the activation or inhibition of a species using normalized Hill functions. The model was used to explore the possible mechanisms for dysregulated signaling in both macrophages and glomerular endothelial cells during diabetic kidney disease progression. The model simulations were trained and validated against in vitro experimental data.

## Authors
Krutika Patidar<sup>a</sup>,  Ashlee N. Ford Versypt<sup>a,b,c</sup>

<sup>a</sup>Department of Chemical and Biological Engineering, University at Buffalo, The State University of New York, Buffalo, NY, USA<br/>
<sup>b</sup>Department of Biomedical Engineering, University at Buffalo, The State University of New York, Buffalo, NY, USA<br/>
<sup>c</sup>Institute for Artificial Intelligence and Data Science, University at Buffalo, The State University of New York, Buffalo, NY, USA<br/>

## Manuscript
K. Patidar and A. N. Ford Versypt, Logic-Based Modeling of Inflammatory Macrophage Crosstalk with Glomerular Endothelial Cells in Diabetic Kidney Disease, bioRxiv preprint, 2023. DOI: 10.1101/2023.04.04.535594 [Preprint](https://biorxiv.org/cgi/content/short/2023.04.04.535594)

## Scripts

* call_ODE_model.m This file calls the following code scripts to perform necessary functions.
* networkODE.m This file contains the ODE equations and utility functions (normalized Hill function) and adds parameter constraints as needed.
* networkODE_opt_loadParams.m This file contains a dictionary of default and optimized parameter values of reaction parameters (W, n, EC50), species parameters (y0, ymax, tau), and species names.
* networkODE_run.m This file runs the ODE model and provides plots that show trained and validated predictions against experimental data.
* networkODE_multirun.m This file uses fitted parameters to run Monte Carlo technique, calculates credible intervals, and visualizes the fitted/validated plots and histograms.
* networkODE_error.m This file computes the sum of squared error (SSE) and weighted SSE between model predictions and data.
* sens.m This file performs a local sensitivity analysis on the time constant (tau), reaction weight (W), Hill coefficient (n), and Half effect (EC50) parameter.
* multistart_param_opt.m This file performs a multi-start parameter estimation using a nonlinear optimizer to estimate values for the sensitive parameters in the model. This file also scales the parameters, samples parameter subsets in a given range using Latin hypercube sampling, and provides the standard deviation of the estimates in each run.
* post_sens.m This file performs local sensitivity analysis on the validated model to identify influential species and interactions in the network for further analyses.
* LHS_Call.m This supporting function uses Latin hypercube sampling to create sample subsets of parameters within a given range.
* jbfill.m This supporting function will fill a region with a color between the two vectors.

## Data
* invitro_data.mat This data file provides an ordered list of normalized data points from published in vitro experiments. The mat file must be loaded to plot model predictions against data.
* The normalized training and validation data is also provided in MIDAS-formatted CSV files.

## STRIKEGOLDD Scripts
* STRIKE_GOLDD.m and options.m These files are provided as part of the STRIKE-GOLDD package. These two scripts are obtained from [STRIKE-GOLDD GitHub Repo](https://github.com/afvillaverde/strike-goldd). These files are used to run the structure identifiability and observability analysis.
* z_network_full_k.m, z_network_full_n.m, z_network_full_W.m, and z_network_full_tau.m These files are used to initialize the network model to perform identifiability analysis for unknown parameters k, n, tau, and W. The structure of these scripts is based on STRIKE-GOLDD package requirements.  
* Additional information on STRIKE-GOLDD package can be found [here](https://github.com/afvillaverde/strike-goldd)

## UQLab Scripts
* UQLab_network_initialize.m This file is used to initialize the network model as per UQLab sensitivity analysis package requirements.
* UQLab_network_run.m This file can be used to run the Sobol sensitivity analysis, provided UQLab package and supporting files are downloaded.
* UQLab_Sobolplots_networkmodel.m This file can be used to plot the First-order and Total order Sobol sensitivity indices for W, n, EC50, and tau parameters.
* Additional information on downloading and using UQLab can be found [here](https://www.uqlab.com/download)

## Recommended Supplementary Packages
* [Netflux](https://github.com/saucermanlab/Netflux) is a package that generates equations and utility functions for networkODE.m

## Acknowledgements
Research reported in this publication was supported by the National Institute of General Medical Sciences of the National Institutes of Health under award number R35GM133763 and NSF CAREER
2133411. The content is solely the responsibility of the authors and does not necessarily represent the official views of the funding agencies.
