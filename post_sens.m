% Author: Krutika Patidar
% Description: Sensitivity analysis of validated model

function [s_FD_Ym, s_FD_W, s_FD_k] = post_sens(params, y0, tspan, sens_change)
tic;
% parameter size Initialization
RP = length(params{1}(1,:));
SP = length(params{3}(:));

% sensitivity coefficient initialization
% size: length(species), length(params), length(time)
s_FD_Ym = zeros(30,SP,1); % Sensitivity index
s_FD_W = zeros(30,RP,1);
s_FD_k = zeros(30,RP,1);

sens_change = 15;
percent = sens_change; % percent perturbation
[t, y] = networkODE_run(tspan,y0,params,0);

% y = y(:,[12, 20, 23, 24, 25, 27, 29, 30]);
% AUC_y = trapz(t, y(:,[12, 18, 20, 23, 25, 27, 29, 30]));


% Parameter: Ymax
for m = 1:SP
    dp = params;
    dp{3}(m) = dp{3}(m)*(1-percent*1e-2); % perturb each "Ymax" parameter by a small amount
    
    [time,dy_model] = networkODE_run(tspan,y0,dp,0); % simulate with perturbed value
    dy_model = real(dy_model);
    dy_modelR = [dy_model(:,[1:30])];
    % AUC_dy = trapz(time, dy_model(:,[12, 13, 20, 23, 24, 25, 27, 29]));
    for l = 1:size(dy_modelR,2)
        %s_FD_tau(l,m,:) = (AUC_dy(l) - AUC_y(l))/params{3}(m)/(percent*1e-2); % difference in AUC values 
        s_FD_Ym(l,m,:) = (dy_modelR(end,l) - y(end,l))/params{3}(m)/(percent*1e-2); % difference in values at 48 hours
    end
end
% Parameter: W
for m = 1:RP
    dp = params;
    dp{1}(1,m) = dp{1}(1,m)*(1-percent*1e-2); % perturb each "W" parameter by a small amount
    
    [time,dy_model] = networkODE_run(tspan,y0,dp,0); % simulate with perturbed value
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
    
    [time,dy_model] = networkODE_run(tspan,y0,dp,0); % simulate with perturbed value
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
%% HEATMAP
%name2 = {'ROSec', 'NO', 'IL-6', 'TNF-\alpha', 'IL-1\beta', 'VEGF-A', 'Calcium', 'Gap Width'};
reac_names = {'=> GLU','=> LPS','LPS => TLR4','GLU => AGE', 'AGE => RAGE','RAGE => NADPH','TLR4 & ROS => NF\kappaB', 'TLR4 => PI3K','NADPH => ROS', 'PI3K => AKT', 'PI3K => ROS', 'NF\kappaB_e_c => TNF-\alpha','AKT => NF\kappaB','NF\kappaB => IL-6','NF\kappaB => TNF-\alpha','NF\kappaB => VEGF-A_m_R_N_A', 'VEGF-A_m_R_N_A => VEGF-A','NF\kappaB => IL-1\beta','VEGF-A => VEGFR1','VEGF-A => VEGFR2','AGE => RAGE_e_c','RAGE_e_c => NADPH_e_c','VEGFR2 => PI3K_e_c','VEGFR1 => PI3K_e_c', 'NADPH_e_c => ROS_e_c', 'PI3K_e_c => AKT_e_c', 'AKT_e_c => eNOS' , 'VEGFR1 => PLC-\gamma' ,'PLC-\gamma => NF\kappaB_e_c' , 'ROS_e_c => NF\kappaB_e_c','NF\kappaB_e_c => IL-6', 'NF\kappaB_e_c => IL-1\beta', 'eNOS  => NO', 'eNOS => ROS_e_c' ,'ROS_e_c & NO => ONOO','!NO => Ca','PLC-\gamma => Ca', 'Ca => pJunction','pJunction => Gap Width', 'Ca => NO'};

figure(80); heatmap(s_FD_Ym(:,:,1), 'Colormap', jet); ax = gca; ax.YData = params{4}(:); ax.XData = params{4}(:);
figure(90); heatmap(s_FD_W(:,:,1), 'Colormap', jet); ax = gca; ax.YData = params{4}(:); ax.XData = reac_names(:);
figure(100); heatmap(s_FD_k(:,:,1), 'Colormap', jet); ax = gca; ax.YData = params{4}(:); ax.XData =  reac_names(:);

%% Species Dyanmics upon perturbation in selected sensitive params
reac_names = {'=> GLU','=> LPS','LPS => TLR4','GLU => AGE', 'AGE => RAGE','RAGE => NADPH','TLR4 & ROS => NFKB', 'TLR4 => PI3K','NADPH => ROS', 'PI3K => AKT', 'PI3K => ROS', 'NFKBec => TNFa','AKT => NFKB','NFKB => IL6','NFKB => TNFa','NFKB => VEGFamRNA', 'VEGFamRNA => VEGFa','NFKB => IL1b','VEGFa => VEGFRec1','VEGFa => VEGFRec2','AGE => RAGEec','RAGEec => NADPHec','VEGFRec2 => PI3Kec','VEGFRec1 => PI3Kec', 'NADPHec => ROSec', 'PI3Kec => AKTec', 'AKTec => eNOS' , 'VEGFRec1 => PLC' ,'PLC => NFKBec' , 'ROSec => NFKBec','NFKBec => IL6', 'NFKBec => IL1b', 'eNOS  => NO', 'eNOS => ROSec' ,'ROSec & NO => ONOO','!NO => Calcium','PLC => Calcium', 'Calcium => pJunc','pJunc => GapWidth', 'Calcium => NO'};

W_pert = [4,16,17,20,21,22,28]; %[26,28,31,33,37,38]; % [4,16,17,20,21,22,24]; %
Ym_pert = [1, 3, 4, 7, 11, 12, 19, 20, 22, 26, 28, 29]; 
Y_spec = [4, 5, 12, 18, 20, 21, 22, 26,28, 29, 30];

figure(5)
for j = 1:length(Ym_pert)
            dp = params;
            %dp{1}(1,W_pert(j)) = dp{1}(1,W_pert(j))*(1 - 50*1e-2); % perturb each parameter by a small amount
            dp{3}(Ym_pert(j)) = dp{3}(Ym_pert(j))*(1 - 50*1e-2);

            options = [];
            [dt, dy] = ode23s(@networkODE,tspan,y0,options,dp);
            for i = 1:11
                hold on
                subplot(4,3,i)
                plot(dt, dy(:,Y_spec(i)), '--', 'LineWidth', 2), xlabel('Time(hour)'),  ylabel(params{4}(Y_spec(i)))
                set(gca,'FontName','Arial','FontSize',12);
                grid on;
                hold on
            end
            

            %legendInfo{j} = num2str(Y_pert(j));

            

                        
end
hold on
%legend(legendInfo)
%legend(reac_names(W_pert))
legend(params{4}(Ym_pert))
end




    
