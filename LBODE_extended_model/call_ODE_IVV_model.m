% call_ODE_IVV_model script to perform the below listed tasks.
% Author: Krutika Patidar
% Description: This file calls various functions ...

%% User-defined input%
clear all

w = warning ('off','all');
% Global data variables
N_ctrl = load('data/finch_density_ctrl.mat');           % Fenestration density for control mice Finch et al. 2022
F_diab = load('data/finch_fenestration_disease.mat');   % Fenestration width/density from Finch et al. 2022
G = load('data/GLU_data.mat');                          % Glucose concentration from Finch et al. 2022 (12-20 weeks) and Lee et al. 2018 (2-11 weeks). Fig 1b male ob-/ob- from both sources

global number_ctrl time_ctrl density diameter GC_conc GC_time GC_LB GC_UB time_lee glu_UB glu_LB time_finch glu_finch LB_lee UB_lee glucose_lee

number_ctrl = N_ctrl.number_ctrl;
time_ctrl = N_ctrl.time_ctrl;
density = F_diab.density;
diameter = F_diab.diameter;
GC_conc = G.GC_conc;
GC_time = G.GC_time;
GC_LB = G.GC_LB;
GC_UB = G.GC_UB;
time_lee = G.time_lee;
glu_UB = G.glu_UB;
glu_LB = G.glu_LB;
time_finch = G.time_finch;
glu_finch = G.glu_finch;
LB_lee = G.LB_lee;
UB_lee = G.UB_lee;
glucose_lee = G.glucose_lee;


state = "diab_mice"; % "diab_mice" for diabetic mice, "norm_mice" for WT mice
prompt = "Choose an option from {""plot_step"", ""MultiStart_NLS"", ""Variability"", ""Predictions"", ""Publication_plots""}: ";
step = input(prompt, 's');

% To run multi-start non-linear least sq optimization, select
% MultiStart_NLS as prompt with necessary modifications for error
% calculation and input parameters.


%%

load data/CGM_db.mat
load data/GLU_data.mat

% Initialize parameters, initial value, and simulation time

% In Vivo Model Params
[params, y0] = coupledODE_IVIVC_params(); %IVV parameters kept separate
[p_params] = coupledODE_physParams(state);


% Time (long-term mice sim.)
    start_time = 2; %weeks
    start_time_h = start_time*7*24;
    end_time = 20; %weeks
    end_time_h = end_time*7*24;
    tspan = start_time_h:1:end_time_h; % hours
% 8000 hours ~ 12 month ~ 48 weeks 

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


%%
if strcmp(step,"plot_step")
    mode = 1; % plots, mode = 0 only simulates
    % Plots from Fitted Model
    if mode>0
        if state == 'norm_mice'
            
            y0(31) = mean(ctrl_finch);
            rng("twister") % Default random number generator algorithm with seed = 0 to ensure that we generate the same sequence of draws
            glu_sampled = zeros(11,1);
            glucose_data_Lee_sd = abs(glucose_lee - LB_lee);
            glucose_data_Finch_sd = abs(glu_finch - glu_UB); 
            for i = 1:length(time_finch)
                glu_sampled(i) = unifrnd(ctrl_LB(:,i), ctrl_UB(:,i)); % 
            end
            [Time, Ypred] = coupledODE_IVV_run(tspan, y0, params, p_params, mode, state, glu_sampled);

        end
        if state == 'diab_mice'
            [Time, Ypred] = coupledODE_IVV_run(tspan, y0, params, p_params, mode, state, GC_conc');
        end
    else
        beep
        'Please provide a correct mode'
    end



elseif strcmp(step,"Publication_plots")
    mode = 1;
    % Plots from Fitted Model
    state = 'diab_mice';
    tau_index = []; W_index = []; n_index = []; k_index = [];
    rng("twister") % Default random number generator algorithm with seed = 0 to ensure that we generate the same sequence of draws
    glu_sampled = zeros(11,1);
    glucose_data_Lee_sd = abs(glucose_lee - LB_lee);
    glucose_data_Finch_sd = abs(glu_finch - glu_UB); 
    subsetIdx = find(time_lee>=6*24*7);
    sigma_data = [glucose_data_Lee_sd(subsetIdx)'  glucose_data_Finch_sd];
    for i = 1:length(GC_time)
        glu_sampled(i) = normrnd(GC_conc(:,i), sigma_data(i)); %
    end
    pub_plots(tspan, y0, params, p_params, mode, state, glu_sampled, tau_index, k_index, n_index, W_index);
    disp('Main-text figures: {2, 3, 4, 5, 6, 7}. Supplementary figures start with 5 and are {52, 53, 54, 55, 56, 57} for B, C, D, E, F, G, respectively.');

elseif strcmp(step, "Predictions")
    disp('1: test_knockout, 2: LSA-based perturbation, 3: time-dependent intervention, 4: Glucose-intervention')
    task = input("choose from {1, 2, 3, 4}: ");
    
    if task < 1 || task > 4
        disp("incorrect input")
    elseif task == 3
        Tstop = input("choose time (t) for intervention (> 336 hr or < 5000 hr): ");
        [s_FD_Ym, s_FD_W] = post_IVV(params, y0, tspan, p_params, state, task, Tstop);
    else
        Tstop = [];
        [s_FD_Ym, s_FD_W] = post_IVV(params, y0, tspan, p_params, state, task, Tstop);
        
    end


elseif strcmp(step, "Variability")
        global_p_best = []; p_fitted = []; error_fitted = [];
        mode = input("choose 3 for glucose variability and 4 for parameter variability "); 
        if mode == 3 || mode == 4
            
            global_p_best = csvread("recal-param\fen_combined_25.csv");
            fitted_p = csvread("recal-param\fen_combined_fitted_25.csv");
            p_fitted = fitted_p(:,1:end-1);
            error_fitted = fitted_p(:,end);
            [Time, Ymean] = coupledODE_IVV_multirun(tspan, y0, params, p_params, mode, state, global_p_best, p_fitted, error_fitted);
        else
            disp("incorrect user input")
        end
         

elseif strcmp(step,"Multistart_NLS")
    % Call Multi-start parameter estimation
    % global_p_best: returns 1xn best parameter set
    % p_fitted: returns 100xn fitted parameter sets
    % error_fitted: returns nx1 fitted sum of squared error
    % *n is the number of parameters for each treatment condition

    % Run parallel specifications
    pLOOP = parpool(8);
    pLOOP.IdleTimeout = 60*8;
    
    repeats = 25;               % integer, number of optimization runs, must be greater than 0
    mode = 0;
    tau_index = []; W_index = []; n_index = []; k_index = [];
    glu_sampled = GC_conc;
    [global_p_best, p_fitted, error_fitted] = multistart_param_opt(params, y0, tspan, p_params, tau_index, n_index, k_index, W_index, repeats, state, glu_sampled);    
    disp(global_p_best);

    delete(pLOOP);

    % Uncomment to calculate minimum RMSE for best set of parameters

    % min_error returns Sum of Sqared Error (SSE)
    [min_error] = coupledODE_IVV_SSE(global_p_best, params, y0, tspan, p_params, state, tau_index, W_index, k_index, glu_sampled);
    disp(min_error);
   
       % p_params(1,2) = global_p_best(1); % EG0
       % p_params(1,9) = global_p_best(2); % Kis
       % p_params(1,12) = global_p_best(3); % Ph
       % p_params(1,14) = global_p_best(4); % w0
       % p_params(1,15) = global_p_best(5); % w
       % p_params(1,16)=  global_p_best(6); % theta
       % p_params(1,17) = global_p_best(7); % At

    if mode==1
        
    [Time, Ypred] = coupledODE_IVV_run(tspan, y0, params, p_params, mode, state, glu_sampled);
    end

end




