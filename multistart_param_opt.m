function [global_p_best, p_fitted, error_fitted] = multistart_param_opt(params, y0, tspan, tau_index, n_index, k_index, W_index)
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
size_W = size(W_index, 2);
%% Scale Parameters (Normalize w.r.t highest value):
%for s = 1: length(params{2})
%    z_tau(s) = (params{2}(s) - min(params{2}))/(max(params{2})-min(params{2})); 
%end
%for s = 1:length(params{1}(2,:))
%    z_n(s) = (params{1}(2,s) - min(params{1}(2,:)))/(max(params{1}(2,:)) -  min(params{1}(2,:)));
%end

z_rpar = [params{1}(1,:); params{1}(2,:); params{1}(3,:)];
z_params = {z_rpar, params{2}(:), params{3}(:), params{4}(:)};                    % Parameters within 0-1 range are not normalized

%% Initialize parameter guess list, create a custom list of parameters to be optimized from the original list of parameters


for tau_iter = 1:size_tau % z_params below
    p_tau(tau_iter) = z_params{2}(tau_index(tau_iter)); % obtained from local_sensitivity analysis
end
for w_iter = 1:size_W  %EC50
    p_W(w_iter) = z_params{1}(1,W_index(w_iter));
end
for n_iter = 1:size_n
    p_n(n_iter) = z_params{1}(2,n_index(n_iter));
end
for k_iter = 1:size_k  %EC50
    p_k(k_iter) = z_params{1}(3,k_index(k_iter));
end


% p_init: parameter guess list

p_init = horzcat(p_tau); %, p_W, p_n, p_k);                        % only valid when atleast 1 sensitive parameter from each category exists


%% Options for parameter estimation with "lsqcurvefit"
fminoptions = optimoptions(@fmincon,'MaxFunctionEvaluations',10000, 'OptimalityTolerance', 1.0000e-12, 'FunctionTolerance', 1.00e-6, 'ConstraintTolerance', 1.000e-12, 'StepTolerance',1e-6);

%fminoptions2 = optimoptions(@fmincon,'MaxFunctionEvaluations',1000, 'OptimalityTolerance', 1.0000e-8, 'FunctionTolerance', 1.00e-8, 'ConstraintTolerance', 1.000e-8);
% global opt solver: MultiStart
% opts = optimoptions(@fmincon,'Algorithm','interior-point');

%% Options for sampling
% total_variables = size(exp_data,2);
total_parameters = size(p_init,2);

% Settings for the multistart optimization
repeats = 100; % more than 1 when necessary

%% Initialize data structures 

% fitted coefficients, Yout for each run, and resnorm (returns the value of the squared 2-norm of the residual at fitted coefficients)
time_for_samples = [0:48];                                                 % Ycalc evaluated at time points
% time_for_samples = [0:28:6272];                                          % Ycalc evaluated at time points 
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
% [Tout, Yout] = ode23(@coupledODE_IVV_hres,tspan,y0,options,params);
% Yout_sampled(1,:,:) = Yout;
p_sampled(1,:) = p_init;
res_norm = networkODE_error(p_init, params, y0, tspan, tau_index, n_index, k_index, W_index) % min_error should be called here


% res_norm = coupledODE_IVV_error(p_init, params, y0, tspan, tau_index, n_index, k_index); 

error_sampled(1) = res_norm


%% Random Sampling

%% Latin Hypercube sampling
% initial parameters  are obtained from local sensitivity analysis
    A = []; B = []; A_eq = []; B_eq = []; % no inequality and equality constraints 
    LB_n = ones(1,size_n); LB_n(1,:) = 1.4; % n>1 required
    UB_n = ones(1, size_n); UB_n(1,:) = 4; % n ~ 3 was most ideal and chosen as default instead of 1.4
    LB_tau = ones(1,size_tau); LB_tau(1,:) = 0.8;% 0.1;
    UB_tau = zeros(1, size_tau); UB_tau(1,:) = 1.2; %10; 
    LB_k = ones(1, size_k); LB_k(1,:) = 0.001; % EC50 0.01?
    UB_k = ones(1,size_k); UB_k(1,:) = 0.84; % EC50 (0.84 calculated from max value of n=4, such that EC50 < 2^-1/n)
    LB_w = zeros(1,size_W); LB_w(1,:) = 0.9; UB_w = ones(1,size_W); %(0,75,1)?
   
    %UB_k(1,5) = 0.01;
    a = horzcat(LB_tau);%, LB_w, LB_n, LB_k);
    b = horzcat(UB_tau);%, UB_w, UB_n, UB_k);


for i = 1:length(p_init)
    p_sampled(2:end,i) = LHS_Call(a(i), p_init(i), b(i), 0 ,repeats,'unif'); % baseline = 10
end

p_guess_rand = p_sampled(2:end,:);
SampledStartPoints = CustomStartPointSet(p_sampled);

%% Check for feasible sampled parameter subsets
 %check_params = params;
 %for i = 2:length(p_sampled(:,1))
 %for m = 1:size_tau
 %       check_params{2}(tau_index(m)) = p_sampled(i,m); %*(max(params{2})-min(params{2})) + min(params{2});        
 %end
 %for n = 1:size_n
 %       check_params{1}(2,n_index(n)) = p_sampled(i, n + size_tau)*(max(params{1}(2,:)) -  min(params{1}(2,:))) + min(params{1}(2,:));
 %end
 %for o = 1:size_k %EC50
 %       check_params{1}(3,k_index(o)) = p_sampled(i, o + size_tau + size_n);
 %end



 %Run coupledODE function with optimal parameter set, no plots or plots training or
 %validation plots based on choice mode = 0, 1, or 2
 %disp(i);
 %[Time, Y_pred] = coupledODE_run(tspan, y0, params, 0);
 
 %options = [];
 
 %[Time, Yout] = ode23s(@coupledODE_IVV_h,tspan,y0,options,params);

 %figure(i)
 %plot(Time, Yout(:,[23, 24, 25, 20,27]));
 %legend('IL-6', 'TNF-a', 'IL-1b', 'NO', 'VEGF')

 %end

%%
% start parallel pool
pLOOP = parpool(8);
pLOOP.IdleTimeout = 60*8;
 
%% Parameter estimation module hereon

tic
parfor rep = 1:repeats % run each of these in parallel
    
    %disp(rep);
    
    F_errmin = @(p) networkODE_error(p, params, y0, tspan, tau_index, n_index, k_index, W_index);
    
    %if isreal(F_errmin(p_sampled(rep+1,:))) == 0
    %    continue
    %end
    
    % fmincon optimization   
    [p_best(rep,:),Obj_best] = fmincon(F_errmin, p_sampled(rep+1,:), A,B,A_eq,B_eq,a,b,[],fminoptions);
    
    p_fitted(rep,:) = p_best(rep,:);
    error_fitted(rep,:) = Obj_best; 
   
    fprintf('run %i finished\n', rep)
    
end
toc

 % MultiStart Problem
    %F_errmin = @(p) networkODE_error(p, params, y0, tspan, tau_index, n_index, k_index, W_index);
    
    
    %opts = fminoptions;
    %problem = createOptimProblem('fmincon','objective',F_errmin,'x0', p_sampled(2,:),'lb',a,'ub',b,'options',opts);
    %ms = MultiStart('UseParallel','always','Display','iter');
    %[p_best,Obj_best] = run(ms,problem,SampledStartPoints);
    %p_fitted = p_best;
    %error_fitted = Obj_best;
    
    %[t, y] = coupledODE_run(tspan, y0, params, 0);
    % res_norm = coupledODE_error(p_best(k,:), params, y0, tspan, tau_index, n_index, k_index); % min_error should be called here
    %Yout_fitted(k,:,:) = y;
 
    




[global_Obj_best,index_best] = min(error_fitted);
global_p_best = p_fitted(index_best,:);

delete(pLOOP);


%writematrix([[W_index]; global_p_best],'params_global/W_best2.csv');
%writematrix([W_index, "error_fitted"; p_fitted, error_fitted],'params_global/W_fitted2.csv');

%writematrix([[n_index]; global_p_best],'params_global/n_best.csv');
%writematrix([n_index, "error_fitted"; p_fitted, error_fitted],'params_global/n_fitted.csv');

%writematrix([[k_index]; global_p_best],'params_global/k_best.csv');
%writematrix([k_index, "error_fitted"; p_fitted, error_fitted],'params_global/k_fitted.csv');

%writematrix([[tau_index]; global_p_best],'params_global/tau_best.csv');
%writematrix([tau_index, "error_fitted"; p_fitted, error_fitted],'params_global/tau_fitted.csv');

writematrix([tau_index, W_index, n_index, k_index; global_p_best],'params_global/LPS_best_tau3.csv');
writematrix([tau_index, W_index, n_index, k_index, "error_fitted"; p_fitted, error_fitted],'params_global/LPS_fitted_tau3.csv');

Params_error = [error_fitted p_fitted]; % coef_fitted is the parameter output from all the LHS fmincon runs
error_column = 1; % first column stores error_fitted
[sorted_Params_error, index] = sortrows(Params_error,error_column);

%-- generate figure to show the error of all the parameter values from LHS
for i = 1:length(p_init)
    plot(1:repeats,sorted_Params_error(:,error_column),'o')
    hold on
end

% savefig("params_global/LPS_minerror_T2.fig")


end

%matlabpool close
