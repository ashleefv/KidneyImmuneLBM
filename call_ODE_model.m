% call_ODE_model script to perform the below listed tasks.
% Author: Krutika Patidar
% Description: This call file calls various functions to run networkODE
% immune response model between macrophages and endothelial cells,
% estimates unknown parameters, optimizes model against experimental
% data from multiple sources, and plots the response variables. 

% Most supporting file names for this model begin with networkODE_###.m
%%
choice = "both"; % options: "GLU", "LPS", "both" 
% choose choice of params (tau, y0, ymax parameters) list based on what combination of stimulus is ON 

% Initialize parameters, initial value, and simulation time
% In Vitro Model Parameters
[params, y0] = networkODE_opt_loadParams(choice); %in vitro as optimized

% Time (short-term sim.)
tspan = [0:1:48];  % Time in hours

sens_change = 0.01; % change in percent
params{1}(1,1) = 1; % change whether GLUCOSE present or absent (either 0 or 1)
params{1}(1,2) = 1; % change whether LPS present or absent (either 0 or 1)


% provide array of indices of tau, n, k parameters that need to be optimized
tau_index = []; n_index = []; k_index = [];

%% Plots from Trained and Valided Model
% calculate confidence intervals within '#_run.m' function

% @kde function documentation: https://www.ics.uci.edu/~ihler/code/kde.html
% @kde function evaluates kde estimates to compute confidence region
% Download kdefunc to use the below functionality

% addpath("~/kdefunc") 

% Run networkODE function with optimal parameter set
% choose mode=0 :no plots
% mode = 1 :plots
% mode = 2 :additional validation plots
% mode = 3 :regulatory response plots
load invitro_data.mat
mode = 2;
[Time, Y_pred] = networkODE_run(tspan, y0, params, mode);
%%
%% Model Workflow
%% Call local sensitivity analysis
% s_FD_tau returns sensitivity coefficients for tau parameter
% s_FD_n returns sensitivity coefficients for n parameter 
% s_FD_k returns sensitivity coefficients for k parameter (EC_50)
% tau_index: index of tau parameter (corresponds to species ID) to be optimized
% n_index: index of n parameter (corresponds to reaction index) to be optimized
% k_index: index of EC_50 parameter (corresponds to reaction index) to be optimized
% params: set of parameters (W, n, EC50, tau, ymax, speciesnames)
% y0: Initial value
% tspan: simulation time points

[s_FD_tau, s_FD_n, s_FD_k, tau_index, n_index, k_index] = sens(params, y0, tspan, sens_change);

disp(tau_index);
disp(n_index);
disp(k_index);
%% Call parameter estimation function

% global_p_best returns set of best parameters
% params: set of parameters (W, n, EC50, tau, ymax, speciesnames)
% y0: Initial values
% tspan: simulation time points
% tau_index: index of tau parameters to be optimized
% n_index: index of n parameters to be optimized
% k_index: index of EC_50 parameter (corresponds to reaction index) to be optimized
% y0_index: index of y0

[global_p_best, h, global_SD] = multistart_param_opt(params, y0, tspan, tau_index, n_index, k_index, y0_index);

disp(global_p_best);
%% Calculate minimum RMSE for best set of parameters

% min_error returns Sum of Sqared Error (SSE)
% w_error returns weighted SSE
% global_p_best returns set of best parameters
% params: set of parameters (W, n, EC50, tau, ymax, speciesnames)
% y0: Initial values
% tspan: simulation time points
% tau_index: index of tau parameters to be tuned
% n_index: index of n parameters to be tuned
% k_index: index of EC_50 parameter (corresponds to reaction index) to be tuned
% y0_index: index of y0

[min_error, w_error] = networkODE_error(global_p_best, params, y0, tspan, tau_index, n_index, k_index);

%% Plot graph of select response variables
% Time: Time output
% Y_pred: Predicted output
% global_p_best: partial list of parameters tuned by fmincon
% params: Complete list of parameters required to run networkODE model
% tau_index, n_index: indices for tau and n parameter list that are tuned
% which are replaced in the complete list of parameters (params)
% y0: Initial value/state
% tspan: simulation time
% mode = 0, does not show plot
% mode = 1, plots training graph and validation data
% mode = 2, plots additional validation graph
% mode = 3, plots regulatory nodes in the network
size_tau = size(tau_index,2);
size_n = size(n_index, 2);
size_k = size(k_index,2);

% update the params list with optimal parameter
%for m = 1:size_tau
  %      params{2}(tau_index(m)) = global_p_best(m);       
%end
%for n = 1:size_n
%       params{1}(2,n_index(n)) = global_p_best(n+size_tau);       
%end
%for o = 1:size_k
%        params{1}(3,k_index(o)) = global_p_best(o+size_tau+size_n);       
%end

% scaled params unscaled to actual values.
for m = 1:size_tau
        params{2}(tau_index(m)) = global_p_best(m)*(max(params{2})-min(params{2})) + min(params{2}); 
        % tau_g(m) = global_p_best(m)*(max(params{2})-min(params{2})) + min(params{2}); 
end

for n = 1:size_n
        params{1}(2,n_index(n)) = global_p_best(n + size_tau)*(max(params{1}(2,:)) -  min(params{1}(2,:))) + min(params{1}(2,:));
        % n_g(n) = global_p_best(n + size_tau)*(max(params{1}(2,:)) -  min(params{1}(2,:))) + min(params{1}(2,:));
end

for o = 1:size_k %EC50
        params{1}(3,k_index(o)) = global_p_best(o+size_tau + size_n);
        % k_g(o) = global_p_best(o+size_tau + size_n);
end


% calculate confidence intervals within '#_run.m' function

%uncomment line below to add @kde function to path
% @kde function documentation: https://www.ics.uci.edu/~ihler/code/kde.html
% @kde function evaluates kde estimates to compute confidence region
% addpath("C:/Research/kdefunc") 
addpath("C:/Users/krutikap/Research/kdefunc") 

% Run networkODE function with optimal parameter set
% choose mode=0 :no plots
% mode = 1 :plots
% mode = 2 :additional validation plots
% mode = 3 :regulatory response plots
mode=0;
[Time, Y_pred] = networkODE_run(tspan, y0, params, mode);


% global_p_best: Optimal parameter set obtained
% var: variance of estimated parameters
% SD: standard deviation or sigma
% OUTPUTS
% CI_lower: Lower bound of CI for each parameter
% CI_upper: Upper bound of CI for each parameter
% var = diag(inv(h));
% SD = sqrt(var);

CI_lower = global_p_best - global_SD*tinv(0.90,length(global_p_best));
CI_upper = global_p_best + global_SD*tinv(0.90,length(global_p_best));





