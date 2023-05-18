% Author: Krutika Patidar
% Description: Sensitivity analysis on three parameters tau, n, k (or
% EC50). Here, available data for response variable IL-6, IL1b, TNFa, VEGFa,
% ROS, eNOS, NO is used to calculate sensitivities by finite difference
% method for each of these parameters. Bar plots can be plotted for each
% variable and sensitivities at last time point.

function [s_FD_tau, s_FD_n, s_FD_k, tau_index, n_index, k_index] = sens(params, y0, tspan, sens_change)
tic;

[t,y] = networkODE_run(tspan,y0,params,0);
yR = y([1,end], [23, 25, 24, 6, 13]);

% Data and parameter size Initialization
c = length(params{1}(1,:));

% tau = params{2};
% n = params{1}(2,:);
% EC50 = params{1}(3,:);


S_FD = zeros(length(yR),c,2); %size: length(species), length(params), length(time)
%S_FD_norm = zeros(length(yR),c,2);

s_FD_tau = zeros(length(yR),length(params{2}(:)),2);
s_FD_n = zeros(length(yR),c,2);
s_FD_k = zeros(length(yR),c,2);

%S_norm_tau = zeros(length(yR),length(params{2}(:)),2);
%S_norm_n = zeros(length(yR),c,2);
%S_norm_k = zeros(length(yR),c,2);

percent = sens_change; %percent parameter perturbation


% Parameter: Tau
for m = 1:length(params{2}(:))
    dp = params;
    dp{2}(m) = dp{2}(m)*(1+percent*1e-2);                   % perturb parameter by a small amount
    
    [time,dy_model] = networkODE_run(tspan,y0,dp,0);
    dy_model = real(dy_model);
    dy_modelR = [dy_model([1,48],[23,25,24,6,13])];
    for l = 1:length(yR)
        s_FD_tau(l,m,:) = (dy_modelR(:,l) - yR(:,l))/params{2}(m)/(percent*1e-2);
        %S_norm_tau(l,m,:) = (dy_modelR(:,1)- yR(:,l))/(percent*1e-2)./yR(:,l);
    end
end

% Parameter: n
for m = 1:length(params{1}(2,:))
    dp = params;
    dp{1}(2,m) = dp{1}(2,m)*(1+percent*1e-2);               % perturb m-th parameter by a small amount
    
    [time,dy_model] = networkODE_run(tspan,y0,dp,0);
    dy_model = real(dy_model);
    dy_modelR = [dy_model([1,48],[23,25,24,6,13])];
    for l = 1:length(yR)
        
        s_FD_n(l,m,:) = (dy_modelR(:,1)- yR(:,l))/params{1}(2,m)/(percent*1e-2);
        %S_norm_n(l,m,:) = (dy_modelR(:,1)- yR(:,l))/(percent*1e-2)./yR(:,l);
    end
end

% Parameter: k
for m = 1:length(params{1}(3,:))
    dp = params;
    dp{1}(3,m) = dp{1}(3,m)*(1+percent*1e-2);                % perturb m-th parameter by a small amount
    [time,dy_model] = networkODE_run(tspan,y0,dp,0);  
    dy_modelR = [dy_model([1,48],[23,25,24,6,13])];
    for l = 1:length(yR)
        s_FD_k(l,m,:) = (dy_modelR(:,1)- yR(:,l))/params{1}(3,m)/(percent*1e-2);
        %S_norm_k(l,m,:) = (dy_modelR(:,1)- yR(:,l))/(percent*1e-2)./yR(:,l);
    
    end
end
toc;
t = toc-tic;
%% Calculate limit/threshold magnitude above which significant parameter indices are considered 

limit_tau = max(max(abs(s_FD_tau(:,:,2)),[],2));
limit_n = max(max(abs(s_FD_n(:,:,2)),[],2));
limit_k = max(max(abs(s_FD_k(:,:,2)),[],2));

% Get the most sensitive parameters (indices)
tau_index = unique([find(abs(s_FD_tau(1,:,2))>0.01*limit_tau), find(abs(s_FD_tau(2,:,2))>0.01*limit_tau)]);
n_index = unique([find(abs(s_FD_n(1,:,2))>0.5*limit_n), find(abs(s_FD_n(2,:,2))>0.5*limit_n)]);
k_index = unique([find(abs(s_FD_k(1,:,2))>0.999*limit_k), find(abs(s_FD_k(2,:,2))>0.999*limit_k)]);

%% BAR PLOT of sensitivities 
top_S_tau = s_FD_tau(:,:,2);
top_S_n = s_FD_n(:,:,2);
top_S_k = s_FD_k(:,:,2);

figure(1)
b = bar(top_S_tau);
name1={'tau';'n';'EC50'};
name2={'IL6', 'IL1b','TNFa', 'VEGFa', 'ROS'};
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


end

    
