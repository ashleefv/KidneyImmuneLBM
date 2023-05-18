function [global_p_best, h, global_SD] = multistart_param_opt(params, y0, tspan, tau_index, n_index, k_index, y0_index)
% Partially adapted from ANFV matlab files 12-29-2021
% Multistart parameter estimation with random or Latin hypercube sampling
% 1. evaluate Yout at the sampled parameters
% 2. use the sampled parameters as input guesses for the parameter estimation to determine new values for the fitted parameters
% 3. evaluate Yout at the fitted parameters
% 4. do runs in parallel
% 5. save output Y_pred and fitted params (p_fitted) 

%% Parameter space initialization
size_tau = size(tau_index,2);
size_n = size(n_index,2);
size_k = size(k_index,2); % EC50
% size_y0 = size(y0_index, 2);

%% Scale Parameters (Normalize w.r.t highest value):
for s = 1: length(params{2})
    z_tau(s) = (params{2}(s) - min(params{2}))/(max(params{2})-min(params{2})); 
end
for s = 1:length(params{1}(2,:))
    z_n(s) = (params{1}(2,s) - min(params{1}(2,:)))/(max(params{1}(2,:)) -  min(params{1}(2,:)));
end

z_rpar = [params{1}(1,:); z_n; params{1}(3,:)];
z_params = {z_rpar, z_tau, params{3}(:), params{4}(:)};                    % Parameters within 0-1 range are not normalized

%% Initialize parameter guess list, create a custom list of parameters to be optimized from the original list of parameters

for m = 1:size_tau % z_params below
    p_tau(m) = z_params{2}(tau_index(m)); % obtained from local_sensitivity analysis
end
for n = 1:size_n
    p_n(n) = z_params{1}(2,n_index(n));
end
for o = 1:size_k  %EC50
    p_k(o) = z_params{1}(3,k_index(o));
end
% for r = 1:size_y0  %EC50
%    p_y0(r) = y0(y0_index(r));
% end

% p_init: parameter guess list

p_init = horzcat(p_tau, p_n, p_k);                                  % only valid when atleast 1 sensitive parameter from each category exists


%% Options for parameter estimation with "lsqcurvefit"
fminoptions = optimoptions(@fmincon,'MaxFunctionEvaluations',3000, 'OptimalityTolerance', 1.0000e-12, 'FunctionTolerance', 1.00e-6, 'ConstraintTolerance', 1.000e-12);

%fminoptions2 = optimoptions(@fmincon,'MaxFunctionEvaluations',1000, 'OptimalityTolerance', 1.0000e-8, 'FunctionTolerance', 1.00e-8, 'ConstraintTolerance', 1.000e-8);
% global opt solver: MultiStart
% opts = optimoptions(@fmincon,'Algorithm','interior-point');

%% Options for sampling
% total_variables = size(exp_data,2);
total_parameters = size(p_init,2);

% Settings for the multistart optimization
repeats = 1; % more than 1 when necessary

%% Initialize data structures 

% fitted coefficients, Yout for each run, and resnorm (returns the value of the squared 2-norm of the residual at fitted coefficients)
time_for_samples = [0:48];                                                 % Ycalc evaluated at time points
length_time = length(time_for_samples); 
p_sampled = zeros(repeats+1, length(p_init));                              % sampled coefficients including nominal
Yout_sampled = zeros(repeats+1,length_time,length(y0));                    % Yout values at time_for_samples for all sampled coefficients
p_fitted = zeros(repeats,length(p_init));                                  % fitted coefficients
error_fitted = zeros(repeats,1); error_fitted(:,1) = 10;                   % resnorm for fitted coefficients at Xdata
error_sampled = zeros(repeats+1,1);                                        % resnorm for sampled coefficients at Xdata 
Yout_fitted = zeros(repeats,length_time,length(y0));                       % Yout values at time_for_samples for all fitted coefficients
%Obj_best = zeros(repeats);                                                % initialize objective function
p_best = zeros(repeats,total_parameters);                                  % intitialize best parameters list
SD_uq = zeros(length(p_init),repeats);                                     % initialize std. deviation for parameter uncertainty
mean_uq = zeros(repeats, length(p_init));                                  % initialize mean for parameter uncertainty

%% Nominal parameter values determined from estimate or previous estimation (CellNOptR)
% p is new parameter list containing parameters to be optimized

%% Single run for initial parameter estimates
[Tout, Yout] = networkODE_run(tspan, y0, params, 0);

options = [];
p_sampled(1,:) = p_init;
res_norm = networkODE_error(p_init, params, y0, tspan, tau_index, n_index, k_index); % min_error should be called here

error_sampled(1) = abs(res_norm);
%% Random Sampling

%% Latin Hypercube sampling
% initial parameters  are obtained from local sensitivity analysis
    A = []; B = []; A_eq = []; B_eq = []; % no inequality and equality constraints 
    LB_n = ones(1,size_n); LB_n(1,:) = 1.7; 
    LB_tau = ones(1,size_tau); LB_tau(1,:) = 0.1;
    UB_n = ones(1, size_n); UB_n(1,:) = 5;
    UB_tau = ones(1, size_tau); UB_tau(1,:) = 10;
    LB_k = ones(1, size_k); LB_k(1,:) = 0.5; %dependent on EC50
    UB_k = ones(1,size_k); UB_k(1,:) = 0.9; %dependent on EC50
    LB_w = [0.1, 0.1]; UB_w = [1.0, 1.0];
    % LB_y0 = zeros(1,size_y0); UB_y0 = ones(1,size_y0);
    a = horzcat(LB_tau, LB_n, LB_k);
    b = horzcat(UB_tau, UB_n, UB_k);

tic
for i = 1:length(p_init)
    p_sampled(2:end,i) = LHS_Call(a(i), p_init(i), b(i), 0 ,repeats,'unif'); 
end
p_guess_rand = p_sampled(2:end,:);
SampledStartPoints = CustomStartPointSet(p_sampled);

%%
p1 = parpool('local', 4);
%% Parameter estimation module hereon
% lsqoptions = optimoptions(@lsqcurvefit,'Algorithm','trust-region-reflective','TolX', 1e-12, 'TolFun', 1e-12);
for k=1:repeats

    disp(k);
    F_errmin = @(p) networkODE_error(p, params, y0, tspan, tau_index, n_index, k_index);
    
    if isreal(F_errmin(p_sampled(k+1,:))) == 0
        continue
    end
    % fmincon  
    [p_best(k,:),Obj_best,exitflag,output,lambda,grad,h] = fmincon(F_errmin, p_sampled(k+1,:), A,B,A_eq,B_eq,LB,UB,[],fminoptions);
   
    
 
    p_fitted(k,:) = p_best(k,:);
    error_fitted(k,:) = Obj_best; 
    % J(k,:) = jacobian;
    SD_uq(:,k) = sqrt(diag(inv(h)));                                       % h: hessian from fmincon output
    mean_uq(k,:) = p_fitted(k,:);
    
end
toc



[global_Obj_best,index_best] = min(error_fitted);
global_p_best = p_best(index_best,:);
global_SD = SD_uq(:,index_best)';
% writematrix([[tau_index, n_index, k_index, y0_index]; global_p_best; global_SD],'params_treatment#_trial#_opt.txt'); 

%%
delete(p1) 
end

%matlabpool close
