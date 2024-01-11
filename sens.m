% Author: Krutika Patidar
% Description: Sensitivity analysis on three parameters tau, n, k (or
% EC50). Here, available data for response variable IL-6, IL1b, TNFa, VEGFa,
% ROS, eNOS, NO is used to calculate sensitivities by finite difference
% method for each of these parameters. Bar plots can be plotted for each
% variable and sensitivities at last time point.

function [s_FD_tau, s_FD_n, s_FD_k, tau_index, n_index, k_index] = sens(params, y0, tspan, sens_change)
tic;
% parameter size Initialization
RP = length(params{1}(1,:));
SP = length(params{3}(:));

% sensitivity coefficient initialization
% size: length(species), length(params), length(time)
s_FD_tau = zeros(30,SP,1); % Sensitivity index
s_FD_W = zeros(30,RP,1);
s_FD_k = zeros(30,RP,1);
s_FD_n = zeros(30,RP,1);


percent = sens_change; % percent perturbation
options =[];
[t, y] = ode23s(@networkODE,tspan,y0,options,params);

% y = y(:,[12, 20, 23, 24, 25, 27, 29, 30]);
% AUC_y = trapz(t, y(:,[12, 18, 20, 23, 25, 27, 29, 30]));


% Parameter: tau
for m = 1:SP
    dp = params;
    dp{2}(m) = dp{2}(m)*(1-percent*1e-2); % perturb each "tau" parameter by a small amount
    
    [time,dy_model] = ode23s(@networkODE,tspan,y0,options,dp); % simulate with perturbed value
    dy_model = real(dy_model);
    dy_modelR = [dy_model(:,[1:30])];
    % AUC_dy = trapz(time, dy_model(:,[12, 13, 20, 23, 24, 25, 27, 29]));
    for l = 1:size(dy_modelR,2)
        %s_FD_tau(l,m,:) = (AUC_dy(l) - AUC_y(l))/params{3}(m)/(percent*1e-2); % difference in AUC values 
        s_FD_tau(l,m,:) = (dy_modelR(end,l) - y(end,l))/params{2}(m)/(percent*1e-2); % difference in values at 48 hours
    end
end
% Parameter: W
for m = 1:RP
    dp = params;
    dp{1}(1,m) = dp{1}(1,m)*(1-percent*1e-2); % perturb each "W" parameter by a small amount
    
    [time,dy_model] = ode23s(@networkODE,tspan,y0,options,dp); % simulate with perturbed value
    dy_model = real(dy_model);
    dy_modelR = [dy_model(:,[1:30])];
    % AUC_dy = trapz(time, dy_model(:,[12, 13, 20, 23, 24, 25, 27, 29]));
    for l = 1:size(dy_modelR,2)
        %s_FD_W(l,m,:) = (AUC_dy(l) - AUC_y(l))/params{1}(1,m)/(percent*1e-2);
        s_FD_W(l,m,:) = (dy_modelR(end,l) - y(end,l))/params{1}(1,m)/(percent*1e-2);
    end
end

% Parameter: k or EC50
for m = 1:RP
    dp = params;
    dp{1}(3,m) = dp{1}(3,m)*(1-percent*1e-2); % perturb each "W" parameter by a small amount
    
    [time,dy_model] = ode23s(@networkODE,tspan,y0,options,dp); % simulate with perturbed value
    dy_model = real(dy_model);
    dy_modelR = [dy_model(:,[1:30])];
    % AUC_dy = trapz(time, dy_model(:,[12, 13, 20, 23, 24, 25, 27, 29]));
    for l = 1:size(dy_modelR,2)
        %s_FD_k(l,m,:) = (AUC_dy(l) - AUC_y(l))/params{1}(3,m)/(percent*1e-2);
        s_FD_k(l,m,:) = (dy_modelR(end,l) - y(end,l))/params{1}(3,m)/(percent*1e-2);
    end
end

% Parameter: n
for m = 1:RP
    dp = params;
    dp{1}(2,m) = dp{1}(2,m)*(1-percent*1e-2); % perturb each "W" parameter by a small amount
    
    [time,dy_model] = ode23s(@networkODE,tspan,y0,options,dp); % simulate with perturbed value
    dy_model = real(dy_model);
    dy_modelR = [dy_model(:,[1:30])];
    % AUC_dy = trapz(time, dy_model(:,[12, 13, 20, 23, 24, 25, 27, 29]));
    for l = 1:size(dy_modelR,2)
        %s_FD_n(l,m,:) = (AUC_dy(l) - AUC_y(l))/params{1}(3,m)/(percent*1e-2);
        s_FD_n(l,m,:) = (dy_modelR(end,l) - y(end,l))/params{1}(2,m)/(percent*1e-2);
    end
end

toc;
t = toc-tic;
%% Calculate limit/threshold magnitude above which significant parameter indices are considered 

limit_tau = max(max(abs(s_FD_tau(:,:,1)),[],2)); % maximum in each parameter set
limit_n = max(max(abs(s_FD_n(:,:,1)),[],2));
limit_k = max(max(abs(s_FD_k(:,:,1)),[],2));
limit_W = max(max(abs(s_FD_W(:,:,1)),[],2));


% WORK ON THIS
% Get the most sensitive parameters (indices)
tau_index = unique([find(abs(s_FD_tau(24,:,1))>0.1*limit_tau), find(abs(s_FD_tau(12,:,1))>0.1*limit_tau)])
n_index = unique([find(abs(s_FD_n(24,:,1))>0.01*limit_n), find(abs(s_FD_n(12,:,1))>0.01*limit_n)])
k_index = unique([find(abs(s_FD_k(24,:,1))>0.01*limit_k), find(abs(s_FD_k(12,:,1))>0.01*limit_k)])
W_index = unique([find(abs(s_FD_W(24,:,1))>0.01*limit_W), find(abs(s_FD_W(12,:,1))>0.01*limit_W)])

%% BAR PLOT of sensitivities 

vars = [23,24,25,27,13,22,20,12];
top_S_tau = s_FD_tau(vars,:,1);
top_S_n = s_FD_n(vars,:,1);
top_S_k = s_FD_k(vars,:,1);
top_S_W = s_FD_W(vars,:,1);

figure(1)
b = bar(top_S_tau);
name1={'tau';'n';'EC50'};
name2=params{4}(vars);
set(gca,'xticklabel',name2);
set(gca,'FontName','Arial','FontSize',10);
grid on;
ylabel('Local Sensitivity of time constant (tau)','FontName','Arial','FontSize',10);


figure(2)
b = bar(top_S_n);
set(gca,'xticklabel',name2);
set(gca,'FontName','Arial','FontSize',10);
grid on;
ylabel('Local Sensitivity of exponent (n)','FontName','Arial','FontSize',10);

figure(3)
b = bar(top_S_k);
set(gca,'xticklabel',name2);
set(gca,'FontName','Arial','FontSize',10);
grid on;
ylabel('Local Sensitivity of EC50 or K parameter','FontName','Arial','FontSize',10);

figure(4)
b = bar(top_S_W);
set(gca,'xticklabel',name2);
set(gca,'FontName','Arial','FontSize',10);
grid on;
ylabel('Local Sensitivity of W parameter','FontName','Arial','FontSize',10);


end

    
