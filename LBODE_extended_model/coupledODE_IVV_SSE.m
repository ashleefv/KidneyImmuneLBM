% Author: Krutika Patidar

% Description: The function evaluates minimum root mean squared error

function [err_db] = coupledODE_IVV_SSE(p, params, y0, tspan, p_params, state, tau_index, W_index, k_index, glu_sampled)

global number_ctrl time_ctrl density diameter GC_conc GC_time GC_LB GC_UB time_lee glu_UB glu_LB time_finch glu_finch LB_lee UB_lee glucose_lee

% uncomment if model-parameters tau, W, n, k are also perturbed.
% size_tau = size(tau_index,2);
% size_W = size(W_index,2);
% size_n = size(n_index,2);
% size_k = size(k_index,2);

% % update the params list with optimal parameter
%  for tau_iter = 1:size_tau
%          params{2}(tau_index(tau_iter)) = p(tau_iter);       
%  end
% for w_iter = 1:size_W 
%         params{1}(1, W_index(w_iter)) = p(w_iter + size_tau);
% end
% for n_iter = 1:size_n 
%         params{1}(2, n_index(n_iter)) = p(n_iter + size_W + size_tau);
% end
%  for k_iter = 1:size_k 
%          params{1}(3,k_index(k_iter)) = p(k_iter + size_W + size_n + size_tau);
%  end



p_params(4) = p(1);
p_params(5) = p(2);
params{2}(34) = p(3);
p_params(6) = p(4);
p_params(7) = p(5);


opts = [];
%opts = odeset('RelTol',1e-20); %,'MaxStep',1e-16);
intv = "none";
[t, y] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,params,p_params, state, glu_sampled, intv);

Yout = real(y);
Tout = t;

% control data for fenestration density

% pred_in = interp1(Tout, Yout(:,37), time_ctrl(:,1)*24*7); error(1,1) = sum((pred_in - number_ctrl(:,1)).^2); pred_in = []; % ctrl data

% diseased data for width and density
pred_in = interp1(Tout, Yout(:,37), time_ctrl(:,1)*24*7); error(1,1) = sum((pred_in - density).^2); pred_in = [];
pred_in = interp1(Tout, Yout(:,38), time_ctrl(:,1)*24*7); error(1,2) = sum((pred_in - diameter).^2); pred_in = [];
  
 
 err_db = sum(error);



end