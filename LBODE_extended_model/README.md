# LBODE Extended Model
Code for Multi-Cellular Network Model Predicts Alterations in Glomerular Endothelial Structure in Diabetic Kidney Disease

[![DOI]()]()

## Overview
This extended logic-based ODE model predicts the effects of long-term exposure to glucose on pro-inflammatory macrophages and glomerular endothelial cells in diabetic mice. (to edit)

## Authors
Krutika Patidar<sup>a</sup>,  Ashlee N. Ford Versypt<sup>a,b,c</sup>

<sup>a</sup>Department of Chemical and Biological Engineering, University at Buffalo, The State University of New York, Buffalo, NY, USA<br/>
<sup>b</sup>Department of Biomedical Engineering, University at Buffalo, The State University of New York, Buffalo, NY, USA<br/>
<sup>c</sup>Institute for Artificial Intelligence and Data Science, University at Buffalo, The State University of New York, Buffalo, NY, USA<br/>

## Manuscript

## Scripts

* call_ODE_model.m 
* networkODE.m 
* networkODE_opt_loadParams.m 
* networkODE_run.m 
* networkODE_multirun.m
* networkODE_error.m 
* multistart_param_opt.m 
* post_sens.m 
* LHS_Call.m This supporting function uses Latin hypercube sampling to create sample subsets of parameters within a given range.
* jbfill.m This supporting function will fill a region with a color between the two vectors.

## Data

## Recommended Supplementary Packages
* [Netflux](https://github.com/saucermanlab/Netflux) is a package that generates equations and utility functions for networkODE.m

## Acknowledgements
Research reported in this publication was supported by the National Institute of General Medical Sciences of the National Institutes of Health under award number R35GM133763 and NSF CAREER
2133411. The content is solely the responsibility of the authors and does not necessarily represent the official views of the funding agencies.
