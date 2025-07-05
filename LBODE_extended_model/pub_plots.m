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
% % Figs 2, B (52)
% [T, Y] = coupledODE_IVV_run(tspan, y0, params, p_params, mode, state, GC_conc');
% 
% % 
% %% Regulatory plots
% 
% global_p_best = []; p_fitted = []; error_fitted = [];
% mode = 3;
% % Figs 4, C (53)
% if mode == 3
%     [T, Y] = coupledODE_IVV_multirun(tspan, y0, params, p_params, mode, state, global_p_best, p_fitted, error_fitted);
% end
% 
% 
% %% Predictions (Fig - 3, 5, 6, 7, D (54), E (55), F (56), G (57))
% 
% % Fig 6 and Table E
% task = 1; Tstop = end_time_h; % placeholder
% [s_FD_Ym, s_FD_W] = post_IVV(params, y0, tspan, p_params, state, task, Tstop);

% % Fig D (53) and E (53)
% task = 2; Tstop = end_time_h; % placeholder
% [s_FD_Ym, s_FD_W] = post_IVV(params, y0, tspan, p_params, state, task, Tstop);
% % 
% task = 3;
% if task == 3
%     Tstop = end_time_h; % Fig G = figure(57) corresponds to 20 weeks
%     [s_FD_Ym, s_FD_W] = post_IVV(params, y0, tspan, p_params, state, task, Tstop);
% 
%     Tstop = 10*7*24; % Fig F = figure(56) corresponds to 10 weeks
%     [s_FD_Ym, s_FD_W] = post_IVV(params, y0, tspan, p_params, state, task, Tstop);
% 
%     Tstop = 8*7*24; % Fig 7 corresponds to 8 weeks
%     [s_FD_Ym, s_FD_W] = post_IVV(params, y0, tspan, p_params, state, task, Tstop);
% end

task = 4;
% Figs 3, 5
if task ==4 
    Tstop = end_time_h; 
    [s_FD_Ym, s_FD_W] = post_IVV(params, y0, tspan, p_params, state, task, Tstop);
end

end
