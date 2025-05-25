% call_ODE_model script to perform the below listed tasks.
% Author: Krutika Patidar
% Description: This file calls various functions to run networkODE
% immune response model between macrophages and endothelial cells,
% identify sensitive parameters, optimizes model parameters using experimental
% data from multiple sources, estimates parameter uncertainty,
% and plots the response variables, and publication plots.

% Most supporting file names for this model begin with networkODE_###.m
%% User-defined input
clear all;
close all;

% Choose treatment condition: "GLU", "LPS", "both" 
choice = "GLU"; 


% Run the model for different steps

% step = 'sim_step'   % Model simulation, No plots
% step = 'plot_step'         % Creates plots
    % mode = 1               % Fitting plots
    % mode = 2               % Validation plots
    % mode = 3               % Regulatory plots
% step = 'sensitivity_step'    % run sensitivity analysis
    % LS = 1                % local sensitivity
    % LS = 0                % user-defined methods
    % sens_change           % float, percent change in parameter
% step = 'pub_plot_step'     % recreates publication plots
% step = 'multistart_opt_step'  % Multi-start optimization, runs in parallel
    % repeats = 100             % integer, number of optimization runs
% step = 'MC_sim_step'          % Monte Carlo Simulation
    % MCmode = 0               % No plots
    % MCmode = 1               % Fitting plots
    % MCmode = 2               % Validation plots

step = 'pub_plot_step';
if strcmp(step, 'pub_plot_step')
    disp('To generate Figures 4, 7, 10, 13, 14, S4A-C (MATLAB Figs 541, 542, 543), S5A (MATLAB Fig 55), S6A (MATLAB Fig 56),')
    disp('set GLU treatment condition for choice on line 15.')
    disp('To generate Figures 5, 8, 11, S5B (MATLAB Fig 55), S6B (MATLAB Fig 56)')
    disp('set LPS treatment condition for choice on line 15.')
    disp('To generate Figures 6, 9, 12, 15, S5C (MATLAB Fig 55), S6C (MATLAB Fig 56), S7A-B (MATLAB Fig 571 and Fig 572),')
    disp('set both treatment condition for choice on line 15.')
end

% choose mode with 'plot_step'
mode = 1;           % specify either 1, 2, or 3

% choose LS, sens_change with 'sensitivity_step'
LS = 1;             % specify either 0 or 1 
sens_change = -0.01; % percent decrease

% choose with 'multistart_opt_step'
repeats = 100;      % must be greater than 0

% choose with 'MC_sim_step'
MCmode = 0;           % specify MCmode ~ 1: fitted plots or 2: validation plots
%%

% Initialize parameters, initial value, and simulation time
% In Vitro Model Parameters
[params, y0] = networkODE_opt_loadParams(choice); %in vitro as optimized

% Time (short-term sim.)
tspan = [0:1:48];  % Time in hours


% params: dictionary of parameters
    % params{1}(1,:): W
    % params{1}(2,:): n
    % params{1}(3,:): k or EC50
    % params{2}(:): tau
    % params{3}(:): ymax
    % params{4}(:): speciesNames
% y0: initial species values
% tspan: simulation time
% tau_index: index of sensitive time constant (tau) parameters    
% W_index: index of sensitive reaction weight parameter (W)
% n_index: index of sensitive Hill coeff. (n) parameter
% k_index: index of sensitive Half effect (EC50) parameter 

if strcmp(choice,"GLU")
    params{1}(1,1) = 1; % change whether GLUCOSE present or absent (either 0 or 1)
    params{1}(1,2) = 0; % change whether LPS present or absent (either 0 or 1)

    % array of parametersindices obtained from global sensitivity analysis (UQLab)
    tau_index = []; tau_index = [ 2     6     9    12    13    14    16    17    19    20    22    23    24    25    27]; 
    W_index = []; W_index = [6     8     9    11    12    14    15    16    17    18    22    25    31    32];
    n_index = []; n_index = [3     4     5     6      7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23    24    25    26    27    30    31    32    33    34    36]; 
    k_index = [];  k_index = [3     5     6     8     9    11    12    14    15    16    17    18    21    22    25    26    27    30    31    32    33  34]; 
    

elseif strcmp(choice,"LPS")
    params{1}(1,1) = 0; % change whether GLUCOSE present or absent (either 0 or 1)
    params{1}(1,2) = 1; % change whether LPS present or absent (either 0 or 1)

    % array of parametersindices obtained from global sensitivity analysis (UQLab)
    tau_index = []; tau_index = [ 2     6     9    12    13    14    16    17    19    20    22    23    24    25    27]; 
     W_index = []; n_index = []; k_index = []; 

elseif strcmp(choice,"both")
    params{1}(1,1) = 1; % change whether GLUCOSE present or absent (either 0 or 1)
    params{1}(1,2) = 1; % change whether LPS present or absent (either 0 or 1)

    % array of parametersindices obtained from global sensitivity analysis (UQLab)
    tau_index = []; tau_index = [ 2     6     9    12    13    14    16    17    19    20    22    23    24    25    27]; 
    W_index = []; n_index = []; k_index = []; 
elseif strcmp(choice, "OFF")
    params{1}(1,1) = 0;
    params{1}(1,2) = 0;
    tau_index = []; W_index = []; n_index = []; k_index = []; 
else
    beep
    'Please provide a valid choice of "GLU", "LPS", "both", or "OFF"'
end

%%
if strcmp(step,"sim_step")
    mode = 0;
    [Time, Y_pred] = networkODE_run(tspan, y0, params, tau_index, k_index, n_index, W_index, mode);
    
elseif strcmp(step,"plot_step")
    % Plots from Fitted Model
    if mode==1
        
       
        [Time, Y_pred] = networkODE_run(tspan, y0, params, tau_index, k_index, n_index, W_index, mode);
    elseif mode==2
        
       
        [Time, Y_pred] = networkODE_run(tspan, y0, params, tau_index, k_index, n_index, W_index, mode);
    elseif mode==3
        
       
        [Time, Y_pred] = networkODE_run(tspan, y0, params, tau_index, k_index, n_index, W_index, mode);
    else
        beep
        'Please provide a valid mode of "1", "2", or "3"'
    end

elseif strcmp(step,"pub_plot_step")
    [Time, Y_pred] = networkODE_pub_plot(tspan, y0, params, tau_index, k_index, n_index, W_index);

    % to circumvent some exportgraphics issue with the tiff files for Figs
    % 13 and 14, we export them here
    if strcmp(choice,"GLU")
        figure(13)
        filename = 'Fig13';
        hold on
        fig = gcf;
        exportgraphics(fig,[filename, 'big.tiff'],'Resolution',1200)

        figure(14)
        filename = 'Fig14';
        hold on
        fig = gcf;
        exportgraphics(fig,[filename, 'big.tiff'],'Resolution',1200)
    end
elseif strcmp(step, "sensitivity_step")

    % LS = 1 (local)
    % LS = 0 (global)
    % s_FD_tau returns sensitivity coefficients for tau parameter
    % s_FD_W returns sensitivity coefficients for W parameter
    % s_FD_n returns sensitivity coefficients for n parameter 
    % s_FD_k returns sensitivity coefficients for k parameter (EC50)
   
if LS == 1
    [s_FD_tau, s_FD_W, s_FD_n, s_FD_k, tau_index, W_index, n_index, k_index] = networkODE_sens(params, y0, tspan, sens_change);
    sprintf('tau_index: %s', num2str(tau_index));
    sprintf('W_index: %s', num2str(W_index))
    sprintf('n_index: %s', num2str(n_index))
    sprintf('k_index (EC50): %s', num2str(k_index))
    
else
    disp("Run UQ Lab global sensitivity analysis scripts")
end


elseif strcmp(step,"multistart_opt_step")
    % Call Multi-start parameter estimation
    % global_p_best: returns 1xn best parameter set
    % p_fitted: returns 100xn fitted parameter sets
    % error_fitted: returns nx1 fitted sum of squared error
    % *n is the number of parameters in each treatment condition

    % specify parallel loop
    pLOOP = parpool(8);
    pLOOP.IdleTimeout = 60*8;
    

    [global_p_best, p_fitted, error_fitted] = multistart_param_opt(params, y0, tspan, tau_index, n_index, k_index, W_index, repeats);    
    disp(global_p_best);

    delete(pLOOP);

    % Uncomment to calculate minimum RMSE for best set of parameters

    % min_error returns Sum of Sqared Error (SSE)
    % [min_error] = networkODE_error(global_p_best, params, y0, tspan, tau_index, n_index, k_index, W_index);
    % disp (min_error);

   

elseif strcmp(step,"MC_sim_step")
    % p_posterior: acceptable set of parameters (1xn)
    % population:  randomly sampled parameter distribution
    % Yp:          posteriors of prediction using Monte Carlo
    % credible:    2.5 and 97.5 quantiles (credible (95%) intervals)

    [p_posterior,population, Yp, credible] = networkODE_multirun(tspan, y0, params, p_fitted, error_fitted, global_p_best, tau_index, k_index, n_index, W_index, MCmode);

else
    beep
    'Please provide a valid mode of "sim_step", "plot_step", "pub_plot_step", "multistart_opt_step", or "MC_sim_step"'

end




