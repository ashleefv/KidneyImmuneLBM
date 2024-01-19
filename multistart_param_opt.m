function [global_p_best, p_fitted, error_fitted] = multistart_param_opt(params, y0, tspan, tau_index, n_index, k_index, W_index, repeats)
% 
% Multistart parameter estimation with Latin hypercube sampling (Partially adapted from ANFV matlab files 12-29-2021)
% 1. Initialize Parameter List
% 2. Initialize Data Structures
% 3. Compute minimum sum of squared error for initial parameter set
% 4. Evaluate output at the initial parameter sets.
% 5. LHS-sampled parameters
% 6. Check for feasible parameters
% 7. Use sampled parameters as initial values for estimation
% 8. Save global-best fit, 100 fitted, and SSE 

%% Parameter space initialization
size_tau = size(tau_index,2);
size_n = size(n_index,2);
size_k = size(k_index,2); % EC50
size_W = size(W_index, 2);
%% Scale Parameters (Normalize w.r.t highest value):

% scale parameters between 0-1 range (tau, n)
% for s = 1: length(params{2})
%     z_tau(s) = (params{2}(s) - min(params{2}))/(max(params{2})-min(params{2})); 
% end
% for s = 1:length(params{1}(2,:))
%     z_n(s) = (params{1}(2,s) - min(params{1}(2,:)))/(max(params{1}(2,:)) -  min(params{1}(2,:)));
% end

z_rpar = [params{1}(1,:); params{1}(2,:); params{1}(3,:)];
z_params = {z_rpar, params{2}(:), params{3}(:), params{4}(:)};                    

%% Initialize parameter guess list, create a custom list of parameters to be optimized from the original list of parameters


for tau_iter = 1:size_tau % tau
    p_tau(tau_iter) = z_params{2}(tau_index(tau_iter)); 
end
for w_iter = 1:size_W     % W
    p_W(w_iter) = z_params{1}(1,W_index(w_iter));
end
for n_iter = 1:size_n     % n
    p_n(n_iter) = z_params{1}(2,n_index(n_iter));
end
for k_iter = 1:size_k     % EC50
    p_k(k_iter) = z_params{1}(3,k_index(k_iter));
end


% p_init: parameter guess list

p_init = horzcat(p_tau, p_W, p_n, p_k);                       
% only valid when atleast 1 sensitive parameter from each category exists

%% Options for parameter estimation with "lsqcurvefit"
fminoptions = optimoptions(@fmincon,'MaxFunctionEvaluations',10000, 'OptimalityTolerance', 1.0000e-12, 'FunctionTolerance', 1.00e-6, 'ConstraintTolerance', 1.000e-12, 'StepTolerance',1e-6);

%fminoptions2 = optimoptions(@fmincon,'MaxFunctionEvaluations',1000, 'OptimalityTolerance', 1.0000e-8, 'FunctionTolerance', 1.00e-8, 'ConstraintTolerance', 1.000e-8);
% global opt solver: MultiStart
% opts = optimoptions(@fmincon,'Algorithm','interior-point');

%% Options for sampling

total_parameters = size(p_init,2);

% Settings for the multistart optimization
% repeats = 100; % more than 1 when necessary

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
% Obj_best = zeros(repeats);                                               % initialize objective function
p_best = zeros(repeats,total_parameters);                                  % intitialize best parameters list
SD_uq = zeros(length(p_init),repeats);                                     % initialize std. deviation for parameter uncertainty
mean_uq = zeros(repeats, length(p_init));                                  % initialize mean for parameter uncertainty


%% Single run for initial parameter estimates
[Tout, Yout] = networkODE_run(tspan, y0, params, tau_index, k_index, n_index, W_index, 0);

options = [];

% Yout_sampled(1,:,:) = Yout;
p_sampled(1,:) = p_init;
res_norm = networkODE_error(p_init, params, y0, tspan, tau_index, n_index, k_index, W_index); 

error_sampled(1) = res_norm;


%% Latin Hypercube sampling
% initial parameters  are obtained from local sensitivity analysis
    A = []; B = []; A_eq = []; B_eq = [];   % no inequality and equality constraints 
    LB_n = ones(1,size_n); LB_n(1,:) = 1.4; % n>1 required
    UB_n = ones(1, size_n); UB_n(1,:) = 4;  % n ~ 3 was most ideal and chosen as default instead of 1.4
    LB_tau = ones(1,size_tau); LB_tau(1,:) = 0.01; 
    UB_tau = zeros(1, size_tau); UB_tau(1,:) = 10; 
    LB_k = ones(1, size_k); LB_k(1,:) = 0.001; % EC50 
    UB_k = ones(1,size_k); UB_k(1,:) = 0.84;   % EC50 (0.84 calculated from max value of n=4, such that EC50 < 2^-1/n)
    LB_w = zeros(1,size_W); LB_w(1,:) = 0.9; UB_w = ones(1,size_W); 
   
 
    a = horzcat(LB_tau, LB_w, LB_n, LB_k);
    b = horzcat(UB_tau, UB_w, UB_n, UB_k);


for i = 1:length(p_init)
    p_sampled(2:end,i) = LHS_Call(a(i), p_init(i), b(i), 0 ,repeats,'unif'); 
end

p_guess_rand = p_sampled(2:end,:);
SampledStartPoints = CustomStartPointSet(p_sampled);

%% Check for feasible sampled parameter subsets
 % check_params = params;
 % for i = 2:length(p_sampled(:,1))
 % for m = 1:size_tau
 %        check_params{2}(tau_index(m)) = p_sampled(i,m);         
 % end
 % for n = 1:size_n
 %        check_params{1}(2,n_index(n)) = p_sampled(i, n + size_tau)
 % end
 % for o = 1:size_k %EC50
 %        check_params{1}(3,k_index(o)) = p_sampled(i, o + size_tau + size_n);
 % end



 % Run coupledODE function with optimal parameter set, no plots or plots training or
 % validation plots based on choice mode = 0, 1, or 2
 % disp(i);
 % [Time, Y_pred] = coupledODE_run(tspan, y0, params, 0);
 
 
 % figure(i)
 % plot(Time, Yout(:,[23, 24, 25, 20,27]));
 % legend('IL-6', 'TNF-a', 'IL-1b', 'NO', 'VEGF')

 %end

%%
% start parallel pool


%% Parameter estimation module hereon

tic
parfor rep = 1:repeats % run each of these in parallel
    
    % disp(rep);
    
    F_errmin = @(p) networkODE_error(p, params, y0, tspan, tau_index, n_index, k_index, W_index);
    
    if isreal(F_errmin(p_sampled(rep+1,:))) == 0
        continue
    end
    
    % fmincon optimization   
    [p_best(rep,:),Obj_best] = fmincon(F_errmin, p_sampled(rep+1,:), A,B,A_eq,B_eq,a,b,[],fminoptions);
    
    p_fitted(rep,:) = p_best(rep,:);
    error_fitted(rep,:) = Obj_best; 
   
    fprintf('run %i finished\n', rep)
    
end
toc  




[global_Obj_best,index_best] = min(error_fitted);
global_p_best = p_fitted(index_best,:);


writematrix([tau_index, W_index, n_index, k_index; global_p_best],'params_global/global_param_best.csv');
writematrix([tau_index, W_index, n_index, k_index, "error_fitted"; p_fitted, error_fitted],'params_global/fitted_params.csv');

Params_error = [error_fitted p_fitted]; % coef_fitted is the parameter output from all the LHS fmincon runs
error_column = 1;                       % first column stores error_fitted
[sorted_Params_error, index] = sortrows(Params_error,error_column);

%-- generate figure to show the error of all the parameter values from LHS
for i = 1:length(p_init)
    plot(1:repeats,sorted_Params_error(:,error_column),'o')
    hold on
end

% savefig("params_global/TR_minError.fig")


end

%matlabpool close
