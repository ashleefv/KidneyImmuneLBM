function pub_plots(tspan, y0, params, p_params, mode, state, glu_sampled, tau_index, k_index, n_index, W_index)

global number_ctrl time_ctrl density diameter GC_conc GC_time GC_LB GC_UB time_lee glu_UB glu_LB time_finch glu_finch LB_lee UB_lee glucose_lee
close all
GLU = params{1}(1,1);

size_tau = size(tau_index,2);
size_n = size(n_index,2);
size_k = size(k_index,2); % EC50
size_W = size(W_index, 2);

% Time (long-term mice sim.)
    start_time = 2; %weeks
    start_time_h = start_time*7*24;
    end_time = 20; %weeks
    end_time_h = end_time*7*24;
    tspan = [start_time_h:1:end_time_h]; % hours

% %% Fitted plots
% mode=1;
% % Figs 2, 4, B (51)
% [T, Y] = coupledODE_IVV_run(tspan, y0, params, p_params, mode, state, GC_conc');


%% Regulatory plots

% global_p_best = []; p_fitted = []; error_fitted = [];
% mode = 3;
% % Fig C (52)
% if mode == 3
%     [T, Y] = coupledODE_IVV_multirun(tspan, y0, params, p_params, mode, state, global_p_best, p_fitted, error_fitted);
% end


%% Predictions (Fig - 3, 5, 6, 7, D (53), E (54), F (55), G (56) )

% Fig 6
task = 1; Tstop = end_time_h; % placeholder
[s_FD_Ym, s_FD_W] = post_IVV(params, y0, tspan, p_params, state, task, Tstop);

% task = 2; Tstop = end_time_h; % placeholder
% [s_FD_Ym, s_FD_W] = post_IVV(params, y0, tspan, p_params, state, task, Tstop);
% 
% task = 3;
% if task == 3
%     Tstop = end_time_h; % figure 3360 corresponds to 20 weeks
%     [s_FD_Ym, s_FD_W] = post_IVV(params, y0, tspan, p_params, state, task, Tstop);
% 
%     Tstop = 10*7*24; % figure 1680 corresponds to 10 weeks
%     [s_FD_Ym, s_FD_W] = post_IVV(params, y0, tspan, p_params, state, task, Tstop);
% 
%     Tstop = 8*7*24; % figure 1344 corresponds to 8 weeks
%     [s_FD_Ym, s_FD_W] = post_IVV(params, y0, tspan, p_params, state, task, Tstop);
% end
% 
% task = 4;
% 
% if task ==4 
%     Tstop = end_time_h; 
%     [s_FD_Ym, s_FD_W] = post_IVV(params, y0, tspan, p_params, state, task, Tstop);
% end

end
