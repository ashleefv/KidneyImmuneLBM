% Author: Krutika Patidar
% Description: Sensitivity analysis of validated model

function [s_FD_Ym, s_FD_W, s_FD_k, W_pert, Ym_pert] = post_sens(params, y0, tspan, sens_change)
tic;

% parameter size initialization
RP = length(params{1}(1,:));
SP = length(params{3}(:));

% sensitivity coefficient initialization
% size: length(species), length(params), length(time)
s_FD_Ym = zeros(length(params{4}(:)),SP,1); % Sensitivity index
s_FD_W = zeros(length(params{4}(:)),RP,1);
s_FD_k = zeros(length(params{4}(:)),RP,1);

% sens_change = 15;
percent = sens_change; % percent perturbation

options=[];
[t, y] = ode23s(@networkODE,tspan,y0,options,params);
y = real(y);

% y = y(:,[12, 20, 23, 24, 25, 27, 29, 30]);
% AUC_y = trapz(t, y(:,[12, 18, 20, 23, 25, 27, 29, 30]));


% Parameter: Ymax
for m = 1:SP
    dp = params;
    dp{3}(m) = dp{3}(m)*(1-percent*1e-2); % perturb each "Ymax" parameter by a small amount
    
    [time,dy_model] = ode23s(@networkODE,tspan,y0,options,dp); % simulate with perturbed value
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

%% HEATMAP
% name2 = {'ROSec', 'NO', 'IL-6', 'TNF-\alpha', 'IL-1\beta', 'VEGF-A', 'Calcium', 'Gap Width'};
% reac_names = {'=> GLU','=> LPS','LPS => TLR4','GLU => AGE', 'AGE => RAGE','RAGE => NADPH','TLR4 & ROS => NFKB', 'TLR4 => PI3K','NADPH => ROS', 'PI3K => AKT', 'PI3K => ROS', 'NFKBec => TNFa','AKT => NFKB','NFKB => IL6','NFKB => TNFa','NFKB => VEGFamRNA', 'VEGFamRNA => VEGFa','NFKB => IL1b','VEGFa => VEGFRec1','VEGFa => VEGFRec2','AGE => RAGEec','RAGEec => NADPHec','VEGFRec2 => PI3Kec','VEGFRec1 => PI3Kec', 'NADPHec => ROSec', 'PI3Kec => AKTec', 'AKTec => eNOS' , 'VEGFRec1 => PLC' ,'PLC => NFKBec' , 'ROSec => NFKBec','NFKBec => IL6', 'NFKBec => IL1b', 'eNOS  => NO', 'eNOS => ROSec' ,'ROSec & NO => ONOO','!NO => Calcium','PLC => Calcium', 'Calcium => pJunc','pJunc => GapWidth', 'Calcium => NO'};
reac_names = {'\Rightarrow GLU','\Rightarrow LPS','LPS \Rightarrow TLR4','GLU \Rightarrow AGE', 'AGE \Rightarrow RAGE','RAGE \Rightarrow NADPH', 'TLR4 & ROS \Rightarrow NF\kappaB', 'TLR4 \Rightarrow PI3K','NADPH \Rightarrow ROS', 'PI3K \Rightarrow AKT', 'PI3K \Rightarrow ROS', 'NF\kappaB_e_c \Rightarrow TNF-\alpha','AKT \Rightarrow NF\kappaB','NF\kappaB \Rightarrow IL-6','NF\kappaB \Rightarrow TNF-\alpha','NF\kappaB \Rightarrow VEGF-A_m_R_N_A', 'VEGF-A_m_R_N_A \Rightarrow VEGF-A','NF\kappaB \Rightarrow IL-1\beta','VEGF-A \Rightarrow VEGFR1','VEGF-A \Rightarrow VEGFR2','AGE \Rightarrow RAGE_e_c','RAGE_e_c \Rightarrow NADPH_e_c','VEGFR2 \Rightarrow PI3K_e_c','VEGFR1 \Rightarrow PI3K_e_c', 'NADPH_e_c \Rightarrow ROS_e_c', 'PI3K_e_c \Rightarrow AKT_e_c', 'AKT_e_c \Rightarrow eNOS' , 'VEGFR1 \Rightarrow PLC-\gamma' ,'PLC-\gamma \Rightarrow NF\kappaB_e_c' , 'ROS_e_c \Rightarrow NF\kappaB_e_c','NF\kappaB_e_c \Rightarrow IL-6', 'NF\kappaB_e_c \Rightarrow IL-1\beta', 'eNOS  \Rightarrow NO', 'eNOS \Rightarrow ROS_e_c' ,'ROS_e_c & NO \Rightarrow ONOO','!NO \Rightarrow Ca','PLC-\gamma \Rightarrow Ca', 'Ca \Rightarrow pJunction','pJunction \Rightarrow Gap Width', 'Ca \Rightarrow NO'};

figure(16); heatmap(s_FD_Ym(:,:,1), 'Colormap', jet); ax = gca; ax.YData = params{4}(:); ax.XData = params{4}(:);
figure(17); heatmap(s_FD_W(:,:,1), 'Colormap', jet); ax = gca; ax.YData = params{4}(:); ax.XData = reac_names(:);
%figure(100); heatmap(s_FD_k(:,:,1), 'Colormap', jet); ax = gca; ax.YData = params{4}(:); ax.XData =  reac_names(:);

% limit is set at maximum of each parameter set
limit_Ym = max(max(real(s_FD_Ym(:,:,1)),[],2));
limit_W = max(max(real(s_FD_W(:,:,1)),[],2));

%% Species Dyanmics upon perturbation in selected sensitive params
reac_names = {'\Rightarrow GLU','\Rightarrow LPS','LPS \Rightarrow TLR4','GLU \Rightarrow AGE', 'AGE \Rightarrow RAGE','RAGE \Rightarrow NADPH', 'TLR4 & ROS \Rightarrow NF\kappaB', 'TLR4 \Rightarrow PI3K','NADPH \Rightarrow ROS', 'PI3K \Rightarrow AKT', 'PI3K \Rightarrow ROS', 'NF\kappaB_e_c \Rightarrow TNF-\alpha','AKT \Rightarrow NF\kappaB','NF\kappaB \Rightarrow IL-6','NF\kappaB \Rightarrow TNF-\alpha','NF\kappaB \Rightarrow VEGF-A_m_R_N_A', 'VEGF-A_m_R_N_A \Rightarrow VEGF-A','NF\kappaB \Rightarrow IL-1\beta','VEGF-A \Rightarrow VEGFR1','VEGF-A \Rightarrow VEGFR2','AGE \Rightarrow RAGE_e_c','RAGE_e_c \Rightarrow NADPH_e_c','VEGFR2 \Rightarrow PI3K_e_c','VEGFR1 \Rightarrow PI3K_e_c', 'NADPH_e_c \Rightarrow ROS_e_c', 'PI3K_e_c \Rightarrow AKT_e_c', 'AKT_e_c \Rightarrow eNOS' , 'VEGFR1 \Rightarrow PLC-\gamma' ,'PLC-\gamma \Rightarrow NF\kappaB_e_c' , 'ROS_e_c \Rightarrow NF\kappaB_e_c','NF\kappaB_e_c \Rightarrow IL-6', 'NF\kappaB_e_c \Rightarrow IL-1\beta', 'eNOS  \Rightarrow NO', 'eNOS \Rightarrow ROS_e_c' ,'ROS_e_c & NO \Rightarrow ONOO','!NO \Rightarrow Ca','PLC-\gamma \Rightarrow Ca', 'Ca \Rightarrow pJunction','pJunction \Rightarrow Gap Width', 'Ca \Rightarrow NO'};

% Most sensitive rules
% W_pert = [3,1,21,28,37,38,39]; 
% Most sensitive species
% Ym_pert = [4,7,17,26,27,28,29]; 

% Indices for rules (W) and species (Ym_pert) that affect Gap Width (except 30: Gap Width)
W_pert = unique(find(abs(s_FD_W(30,:,1))>0.2*limit_W));
Ym_pert = unique(find(abs(s_FD_Ym(30,[1:29],1))>0.2*limit_Ym));

% disp(W_pert)
% disp(Ym_pert)

Y_spec = [30];
linest = ["-", "--", "-.", ":", "--.",'-','--','-.'];

% plot the effect of 50% knockdown of select reaction rules and species on
% Gap Width change
figure(18)
hold on
for j = 1:length(Ym_pert)
            dp = params;
            dp{3}(Ym_pert(j)) = dp{3}(Ym_pert(j))*(1 - 50*1e-2);

            options = [];
            [dt, dy] = ode23s(@networkODE,tspan,y0,options,dp);
            dy = real(dy);

            hold on
            subplot(1,2,2)
            plot(dt, dy(:,Y_spec(1)), linest(j), 'LineWidth', 2), xlabel('Time (hour)'),  ylabel(append('Change in ',params{4}(Y_spec(1))))
            set(gca,'FontName','Arial','FontSize',16);
            ylim([0,1.05])
            yticks([0,0.2,0.4,0.6,0.8,1])
            xlim([0,48])
            xticks([0 12 24 36 48])
            hold on
           
          
            

                        
end
hold on
legend(params{4}(Ym_pert), 'Location', 'SouthEast')

figure(18)
hold on
for j = 1:length(W_pert)
            dp = params;
            dp{1}(1,W_pert(j)) = dp{1}(1,W_pert(j))*(1 - 50*1e-2); % perturb each parameter by a small amount

            options = [];
            [dt, dy] = ode23s(@networkODE,tspan,y0,options,dp);
            dy = real(dy);
            
            hold on
            subplot(1,2,1)
            plot(dt, dy(:,Y_spec(1)), linest(j), 'LineWidth', 2), xlabel('Time (hour)'),  ylabel(append('Change in ',params{4}(Y_spec(1))))
            set(gca,'FontName','Arial','FontSize',16);
            ylim([0,1.05])
            yticks([0,0.2,0.4,0.6,0.8,1])
            xlim([0,48])
            xticks([0 12 24 36 48])
            hold on
            
                        
end
hold on
legend(reac_names(W_pert), 'Location', 'SouthEast')
end



    
