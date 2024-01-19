% Author: Krutika Patidar
% Description: Sensitivity analysis on three parameters tau, W, n, k (or
% EC50). Here, available data for output variables is used to calculate sensitivities 
% by finite difference method for each of these parameters. 
% Bar plots can be plotted for each variable and sensitivities 
% at end time point (48 hours).

function [s_FD_tau, s_FD_W, s_FD_n, s_FD_k, tau_index, W_index, n_index, k_index] = networkODE_sens(params, y0, tspan, sens_change)
tic;
% parameter size Initialization
RP = length(params{1}(1,:));
SP = length(params{3}(:));

% sensitivity coefficient initialization
% size: length(species), length(params), length(time)
s_FD_tau = zeros(length(params{4}(:)),SP,1); % Sensitivity index
s_FD_W = zeros(length(params{4}(:)),RP,1);
s_FD_k = zeros(length(params{4}(:)),RP,1);
s_FD_n = zeros(length(params{4}(:)),RP,1);


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
        %s_FD_tau(l,m,:) = (AUC_dy(l) - AUC_y(l))/params{2}(m)/(percent*1e-2); % difference in AUC values 
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

% Parameter: n
for m = 1:RP
    dp = params;
    dp{1}(2,m) = dp{1}(2,m)*(1-percent*1e-2); % perturb each "n" parameter by a small amount
    
    [time,dy_model] = ode23s(@networkODE,tspan,y0,options,dp); % simulate with perturbed value
    dy_model = real(dy_model);
    dy_modelR = [dy_model(:,[1:30])];
    % AUC_dy = trapz(time, dy_model(:,[12, 13, 20, 23, 24, 25, 27, 29]));
    for l = 1:size(dy_modelR,2)
        %s_FD_n(l,m,:) = (AUC_dy(l) - AUC_y(l))/params{1}(2,m)/(percent*1e-2);
        s_FD_n(l,m,:) = (dy_modelR(end,l) - y(end,l))/params{1}(2,m)/(percent*1e-2);
    end
end

% Parameter: k or EC50
for m = 1:RP
    dp = params;
    dp{1}(3,m) = dp{1}(3,m)*(1-percent*1e-2); % perturb each "EC50" parameter by a small amount
    
    [time,dy_model] = ode23s(@networkODE,tspan,y0,options,dp); % simulate with perturbed value
    dy_model = real(dy_model);
    dy_modelR = [dy_model(:,[1:30])];
    % AUC_dy = trapz(time, dy_model(:,[12, 13, 20, 23, 24, 25, 27, 29]));
    for l = 1:size(dy_modelR,2)
        %s_FD_k(l,m,:) = (AUC_dy(l) - AUC_y(l))/params{1}(3,m)/(percent*1e-2);
        s_FD_k(l,m,:) = (dy_modelR(end,l) - y(end,l))/params{1}(3,m)/(percent*1e-2);
    end
end


toc;
t = toc-tic;
%% Calculate limit/threshold magnitude above which significant parameter indices are considered 

limit_tau = max(max(real(s_FD_tau(:,:,1)),[],2)); % limit is set as the maximum in each parameter set
limit_n = max(max(real(s_FD_n(:,:,1)),[],2));
limit_k = max(max(real(s_FD_k(:,:,1)),[],2));
limit_W = max(max(real(s_FD_W(:,:,1)),[],2));


% Get the most sensitive parameters (indices) that are within 10% of the limit

tau_index = []; W_index=[]; n_index= []; k_index=[];
for i=1:length(params{4})
    tau_Id = find(abs(s_FD_tau(i,:,1))>0.1*limit_tau);
    if isempty(tau_Id)
        
    else
       
       tau_index = [tau_index, tau_Id];
       
    end
end
tau_index = unique(tau_index);

for i=1:length(params{4})
    W_Id = find(abs(s_FD_W(i,:,1))>0.1*limit_W);
    if isempty(W_Id)
        
    else
        W_index = [W_index, W_Id];
       
    end
end
W_index = unique(W_index);

for i=1:length(params{4})
    n_Id = find(abs(s_FD_n(i,:,1))>0.1*limit_n);
    if isempty(n_Id)
        
    else
     n_index = [n_index, n_Id];
       
    end
end
n_index = unique(n_index);

for i=1:length(params{4})
    k_Id = find(abs(s_FD_k(i,:,1))>0.1*limit_k);
    if isempty(k_Id)
        
    else
        k_index = [k_index, k_Id];
       
    end
end
k_index = unique(k_index);

%% BAR PLOT of sensitivities 

vars = [23,24,25,27,13,22,20,12];
top_S_tau = real(s_FD_tau(vars,:,1));
top_S_n = real(s_FD_n(vars,:,1));
top_S_k = real(s_FD_k(vars,:,1));
top_S_W = real(s_FD_W(vars,:,1));

figure(1)
b = bar(top_S_tau);

% name1={'tau';'n';'EC50'};
name = params{4}(vars);

set(gca,'xticklabel',name);
set(gca,'FontName','Arial','FontSize',10);
grid on;
ylabel('Local Sensitivity of \tau','FontName','Arial','FontSize',10);


figure(2)
b = bar(top_S_n);
set(gca,'xticklabel',name);
set(gca,'FontName','Arial','FontSize',10);
grid on;
ylabel('Local Sensitivity of n parameter','FontName','Arial','FontSize',10);

figure(3)
b = bar(top_S_k);
set(gca,'xticklabel',name);
set(gca,'FontName','Arial','FontSize',10);
grid on;
ylabel('Local Sensitivity of EC_{50} parameter','FontName','Arial','FontSize',10);

figure(4)
b = bar(top_S_W);
set(gca,'xticklabel',name);
set(gca,'FontName','Arial','FontSize',10);
grid on;
ylabel('Local Sensitivity of W parameter','FontName','Arial','FontSize',10);


end

    
