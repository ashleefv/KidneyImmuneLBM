# README #

This file guides you through the coupledODE_extended# model files and their function

## Description of coupledODE_extended#.m files ##

coupledODE_extended.m holds model equations, model parameter constraints.
coupldeODE_extended_loadParams.m holds parameter values.
coupledODE_run.m runs coupledODE equations for time (t), initial conditions, and parameter values. 
coupldeODE_error.m computes the SSE error between model and data.
multistart_param_opt.m performs nonlin-leastsq optimization of parameters.
sens.m computes local sensitivities of parameters.


## Call 
call call_coupldedODE_model.m to perform step-wise functionalities as mentioned above