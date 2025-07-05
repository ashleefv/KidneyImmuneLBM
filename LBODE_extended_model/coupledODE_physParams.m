function [p_params] = coupledODE_physParams(state)

if state == "norm_mice" || state == "diab_mice"

% step 1 (formation based parameters only)
% see recal-param\fenC_step1_global_25.csv
    yss = 7; % unitless
    nf = 4.00; % unitless
    kform = 1.01; % 1/hr
%step 2 (other parameters)
% see recal-param\fen_combined_global_25.csv
    ks = 65.9; % nm/hr
    kd = 2.04; % % 1/hr
    kloss = 4.61; % 1/hr
    yss2 = 4.02; % unitless

end

p_params = [yss, nf, kform, ks, kd, kloss, yss2]; 

end