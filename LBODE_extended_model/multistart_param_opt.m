function [global_p_best, p_fitted, error_fitted] = multistart_param_opt(params, y0, tspan, p_params, tau_index, n_index, k_index, W_index, repeats, state, glu_sampled)
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
%    z_tau(s) = (params{2}(s) - min(params{2}(:)))/(max(params{2}(:))-min(params{2}(:)));
%end
%for s = 1:length(params{1}(2,:))
%    z_n(s) = (params{1}(2,s) - min(params{1}(2,:)))/(max(params{1}(2,:)) -  min(params{1}(2,:)));
%end

% z_rpar = [params{1}(1,:); params{1}(2,:); params{1}(3,:)];
% z_params = {z_rpar, params{2}(:), params{3}(:), params{4}(:)};                    % Parameters within 0-1 range are not normalized

%% Initialize parameter guess list, create a custom list of parameters to be optimized from the original list of parameters


% for m = 1:size_tau % z_params below
%     p_tau(m) = z_params{2}(tau_index(m)); % obtained from local_sensitivity analysis
% end
% 
% for s = 1:size_W  %W
%     p_w(s) = z_params{1}(1,W_index(s));
% end
% 
% % for n = 1:size_n
% %     p_n(n) = z_params{1}(2,n_index(n));
% % end
%  for o = 1:size_k  %EC50
%      p_k(o) = z_params{1}(3,k_index(o));
%  end



%%

%R0 = p_params(1,1); %
%EG0 = p_params(1,2); 
%Si = p_params(1,3); 
%Imax = p_params(1,4);
%alpha = p_params(1,5); 
%ki = p_params(1,6); 
%k1 = p_params(1,7); 
%k2 = p_params(1,8); 
%Kis = p_params(1,9); 
%Vh = p_params(1,10);
%Br = p_params(1,11); 
%Ph = p_params(1,12); 
%Ac = p_params(1,13); 
%w0 = p_params(1,14); 
%w = p_params(1,15);
%theta = p_params(1,16); 
%At =  p_params(1,17); 
%HG_max = p_params(1,18); 
%HG_min = p_params(1,19); 
% Ci = p_params(end);

% p_init: parameter guess list

p_init = horzcat(p_params(4), p_params(5), params{2}(33), p_params(6), p_params(7));                  % order tau, W, k, n, ... 
                                                    % only valid when atleast 1 sensitive parameter from each category exists

disp(p_init);
%% Options for parameter estimation with "lsqcurvefit"
fminoptions = optimoptions(@fmincon,'MaxFunctionEvaluations',3000, 'OptimalityTolerance', 1.0000e-12, 'FunctionTolerance', 1.00e-6, 'ConstraintTolerance', 1.000e-12);

fminoptions2 = optimoptions(@fmincon,'MaxFunctionEvaluations',100, 'OptimalityTolerance', 1.0000e-8, 'FunctionTolerance', 1.00e-8, 'ConstraintTolerance', 1.000e-8);
% global opt solver: MultiStart
% opts = optimoptions(@fmincon,'Algorithm','interior-point');

%% Options for sampling
% total_variables = size(exp_data,2);
total_parameters = size(p_init,2);

% Settings for the multistart optimization
% repeats = 1; % more than 1 when necessary

%% Initialize data structures 

% fitted coefficients, Yout for each run, and resnorm (returns the value of the squared 2-norm of the residual at fitted coefficients)
%time_for_samples = [0:48];                                                 % Ycalc evaluated at time points
time_for_samples = tspan;                                          % Ycalc evaluated at time points 
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
%[Tout, Yout] = coupledODE_run(tspan, y0, params, 0);

options= [];
[Tout, Yout] = coupledODE_IVV_run(tspan, y0, params, p_params, 0, state, glu_sampled);
% Yout_sampled(1,:,:) = Yout;
p_sampled(1,:) = p_init;

%
[res_norm] = coupledODE_IVV_SSE(p_init, params, y0, tspan, p_params, state, tau_index, W_index, k_index, glu_sampled);


error_sampled(1) = res_norm


%% Random Sampling

%% Latin Hypercube sampling
% initial parameters  are obtained from local sensitivity analysis
    A = []; B = []; A_eq = []; B_eq = []; % no inequality and equality constraints 
   

%       a = [0.75, 0.25, 0.01]; % pJunction
%       b = [1, 0.84, 0.25];


%       a = [1, 1, 0.5, 0.5, 0.01, 0.01];
%       b = [1000, 1000, 1, 1, 0.84, 0.84]; % Diameter

        a = [45, 1, 400, 1, 3];
        b = [75, 4, 600, 5, 5];

tic
for i = 1:length(p_init)
    p_sampled(2:end,i) = LHS_Call(a(i), p_init(i), b(i), 0 ,repeats,'unif'); % baseline = 10
end
p_guess_rand = p_sampled(2:end,:);
SampledStartPoints = CustomStartPointSet(p_sampled);

disp(p_sampled);
%% Parameter estimation module hereon

parfor rep=1:repeats
    
   
    
    F_errmin = @(p) coupledODE_IVV_SSE(p, params, y0, tspan, p_params, state, tau_index, W_index, k_index, glu_sampled);
    
    if isreal(F_errmin(p_sampled(rep+1,:))) == 0
        disp(isreal(F_errmin(p_sampled(rep+1,:))));
        continue
    end
    % fmincon  
    [p_best(rep,:),Obj_best] = fmincon(F_errmin, p_sampled(rep+1,:), A,B,A_eq,B_eq,a,b,[],fminoptions);
   
   
 
    p_fitted(rep,:) = p_best(rep,:);
    error_fitted(rep,:) = Obj_best; 
    
    fprintf('OPT-RUN %i finished\n', rep)
    
end
toc


%%

[global_Obj_best,index_best] = min(error_fitted);
global_p_best = p_best(index_best,:);

%%
writematrix([global_p_best],'recal-param/fen_combined_25.csv');
writematrix([ p_fitted, error_fitted],'recal-param/fen_combined_fitted_25.csv');

Params_error = [error_fitted p_fitted]; % coef_fitted is the parameter output from all the LHS fmincon runs
error_column = 1; % first column stores error_fitted
[sorted_Params_error, index] = sortrows(Params_error,error_column);

%-- generate figure to show the error of all the parameter values from LHS
for i = 1:length(p_init)
    plot(1:repeats,sorted_Params_error(:,error_column),'o')
    hold on
end
savefig("recal-param/SSE_fen_combined_error.fig")


end

% matlabpool close
