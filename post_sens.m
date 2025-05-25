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

% sens_change = -15;
percent = sens_change; % percent perturbation

options=[];
[t, y] = ode23s(@networkODE,tspan,y0,options,params);
y = real(y);

% y = y(:,[12, 20, 23, 24, 25, 27, 29, 30]);
% AUC_y = trapz(t, y(:,[12, 18, 20, 23, 25, 27, 29, 30]));


% Parameter: Ymax
for m = 1:SP
    dp1 = params;
    dp1{3}(m) = dp1{3}(m)*(1+percent*1e-2); % perturb each "Ymax" parameter by a small amount
    
    [time,dy_model] = ode23s(@networkODE,tspan,y0,options,dp1); % simulate with perturbed value
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
    dp1 = params;
    dp1{1}(1,m) = dp1{1}(1,m)*(1+percent*1e-2); % perturb each "W" parameter by a small amount
    
    [time,dy_model] = ode23s(@networkODE,tspan,y0,options,dp1); % simulate with perturbed value
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
    dp1 = params;
    dp1{1}(3,m) = dp1{1}(3,m)*(1+percent*1e-2); % perturb each "EC50" parameter by a small amount
    
    [time,dy_model] = ode23s(@networkODE,tspan,y0,options,dp1); % simulate with perturbed value
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

%% HEATMAP (FIG S7A(1) and S7B(2) from supplemental)
% name2 = {'ROSec', 'NO', 'IL-6', 'TNF-\alpha', 'IL-1\beta', 'VEGF-A', 'Calcium', 'Gap Width'};
% reac_names = {'=> GLU','=> LPS','LPS => TLR4','GLU => AGE', 'AGE => RAGE','RAGE => NADPH','TLR4 & ROS => NFKB', 'TLR4 => PI3K','NADPH => ROS', 'PI3K => AKT', 'PI3K => ROS', 'NFKBec => TNFa','AKT => NFKB','NFKB => IL6','NFKB => TNFa','NFKB => VEGFamRNA', 'VEGFamRNA => VEGFa','NFKB => IL1b','VEGFa => VEGFRec1','VEGFa => VEGFRec2','AGE => RAGEec','RAGEec => NADPHec','VEGFRec2 => PI3Kec','VEGFRec1 => PI3Kec', 'NADPHec => ROSec', 'PI3Kec => AKTec', 'AKTec => eNOS' , 'VEGFRec1 => PLC' ,'PLC => NFKBec' , 'ROSec => NFKBec','NFKBec => IL6', 'NFKBec => IL1b', 'eNOS  => NO', 'eNOS => ROSec' ,'ROSec & NO => ONOO','!NO => Calcium','PLC => Calcium', 'Calcium => pJunc','pJunc => GapWidth', 'Calcium => NO'};
reac_names = {'\Rightarrow GLU','\Rightarrow LPS','LPS \Rightarrow TLR4','GLU \Rightarrow AGE', 'AGE \Rightarrow RAGE','RAGE \Rightarrow NADPH', 'TLR4 & ROS \Rightarrow NF\kappaB', 'TLR4 \Rightarrow PI3K','NADPH \Rightarrow ROS', 'PI3K \Rightarrow AKT', 'PI3K \Rightarrow ROS', 'NF\kappaB_e_c \Rightarrow TNF-\alpha','AKT \Rightarrow NF\kappaB','NF\kappaB \Rightarrow IL-6','NF\kappaB \Rightarrow TNF-\alpha','NF\kappaB \Rightarrow VEGF-A_m_R_N_A', 'VEGF-A_m_R_N_A \Rightarrow VEGF-A','NF\kappaB \Rightarrow IL-1\beta','VEGF-A \Rightarrow VEGFR1','VEGF-A \Rightarrow VEGFR2','AGE \Rightarrow RAGE_e_c','RAGE_e_c \Rightarrow NADPH_e_c','VEGFR2 \Rightarrow PI3K_e_c','VEGFR1 \Rightarrow PI3K_e_c', 'NADPH_e_c \Rightarrow ROS_e_c', 'PI3K_e_c \Rightarrow AKT_e_c', 'AKT_e_c \Rightarrow eNOS' , 'VEGFR1 \Rightarrow PLC-\gamma' ,'PLC-\gamma \Rightarrow NF\kappaB_e_c' , 'ROS_e_c \Rightarrow NF\kappaB_e_c','NF\kappaB_e_c \Rightarrow IL-6', 'NF\kappaB_e_c \Rightarrow IL-1\beta', 'eNOS  \Rightarrow NO', 'eNOS \Rightarrow ROS_e_c' ,'ROS_e_c & NO \Rightarrow ONOO','!NO \Rightarrow Ca','PLC-\gamma \Rightarrow Ca', 'Ca \Rightarrow pJunction','pJunction \Rightarrow Gap Width', 'Ca \Rightarrow NO'};


figure(572); heatmap(s_FD_Ym([1:29],[1:29],1), 'Colormap', hot, 'CellLabelColor', 'None'); ax = gca; ax.YData = params{4}([1:29]); ax.XData = params{4}([1:29]);
figure(571); heatmap(s_FD_W([1:29],[1:38,40],1), 'Colormap', hot, 'CellLabelColor', 'None'); ax = gca; ax.YData = params{4}([1:29]); ax.XData = reac_names([1:38,40]);
%figure(100); heatmap(s_FD_k(:,:,1), 'Colormap', jet); ax = gca; ax.YData = params{4}(:); ax.XData =  reac_names(:);

% limit is set at maximum of each parameter set
limit_Ym = max(max(real(s_FD_Ym([1:29],[1:29],1)),[],2));
limit_W = max(max(real(s_FD_W([1:29],[1:38,40],1)),[],2));

%% Species Dyanmics upon perturbation in selected sensitive params
reac_names = {'\Rightarrow GLU','\Rightarrow LPS','LPS \Rightarrow TLR4','GLU \Rightarrow AGE', 'AGE \Rightarrow RAGE','RAGE \Rightarrow NADPH', 'TLR4 & ROS \Rightarrow NF\kappaB', 'TLR4 \Rightarrow PI3K','NADPH \Rightarrow ROS', 'PI3K \Rightarrow AKT', 'PI3K \Rightarrow ROS', 'NF\kappaB_e_c \Rightarrow TNF-\alpha','AKT \Rightarrow NF\kappaB','NF\kappaB \Rightarrow IL-6','NF\kappaB \Rightarrow TNF-\alpha','NF\kappaB \Rightarrow VEGF-A_m_R_N_A', 'VEGF-A_m_R_N_A \Rightarrow VEGF-A','NF\kappaB \Rightarrow IL-1\beta','VEGF-A \Rightarrow VEGFR1','VEGF-A \Rightarrow VEGFR2','AGE \Rightarrow RAGE_e_c','RAGE_e_c \Rightarrow NADPH_e_c','VEGFR2 \Rightarrow PI3K_e_c','VEGFR1 \Rightarrow PI3K_e_c', 'NADPH_e_c \Rightarrow ROS_e_c', 'PI3K_e_c \Rightarrow AKT_e_c', 'AKT_e_c \Rightarrow eNOS' , 'VEGFR1 \Rightarrow PLC-\gamma' ,'PLC-\gamma \Rightarrow NF\kappaB_e_c' , 'ROS_e_c \Rightarrow NF\kappaB_e_c','NF\kappaB_e_c \Rightarrow IL-6', 'NF\kappaB_e_c \Rightarrow IL-1\beta', 'eNOS  \Rightarrow NO', 'eNOS \Rightarrow ROS_e_c' ,'ROS_e_c & NO \Rightarrow ONOO','!NO \Rightarrow Ca','PLC-\gamma \Rightarrow Ca', 'Ca \Rightarrow pJunction','pJunction \Rightarrow Gap Width', 'Ca \Rightarrow NO'};

% Most sensitive rules
% W_pert = [3,1,21,28,37,38,39]; 
% Most sensitive species
% Ym_pert = [4,7,17,26,27,28,29]; 

% Indices for rules (W) and species (Ym_pert) that affect Gap Width (except 30: Gap Width)
W_pert = unique(find(abs(s_FD_W(28,[1:38,40],1))>0.1));
Ym_pert = unique(find(abs(s_FD_Ym(28,[1:29],1))>0.1));

% disp(W_pert)
% disp(Ym_pert)

Y_spec = [28];
linest = ["-", "--", "-.", ":", "--.",'-','--','-.'];

% plot the effect of 50% knockdown of select reaction rules and species on
% pJunc change
figure(15)
options = [];
[T, Y] = ode23s(@networkODE,tspan,y0,options,params);
Y = real(Y);
labelstring = {'a', 'b', 'c', 'd'};

hold on

for j = 1:length(Ym_pert)
            dp1 = params;
            dp1{3}(Ym_pert(j)) = dp1{3}(Ym_pert(j))*(1 - 15*1e-2);

            options = [];
            [dt, dy] = ode23s(@networkODE,tspan,y0,options,dp1);
            dy = real(dy);

            hold on
            subplot(2,2,2)
            plot(dt, dy(:,Y_spec(1)), linest(j), 'LineWidth', 2), xlabel('Time (hour)'),  ylabel(append('Normalized ',params{4}(Y_spec(1))))
            set(gca,'FontName','Arial','FontSize',16);
            ylim([0,1.05])
            yticks([0,0.2,0.4,0.6,0.8,1])
            xlim([0,48])
            xticks([0 12 24 36 48])


            dp2 = params;
            dp2{3}(Ym_pert(j)) = dp2{3}(Ym_pert(j))*(1 - 50*1e-2);

            options = [];
            [dt, dy] = ode23s(@networkODE,tspan,y0,options,dp2);
            dy = real(dy);

            hold on
            subplot(2,2,4)
            plot(dt, dy(:,Y_spec(1)), linest(j), 'LineWidth', 2), xlabel('Time (hour)'),  ylabel(append('Normalized ',params{4}(Y_spec(1))))
            set(gca,'FontName','Arial','FontSize',16);
            ylim([0,1.05])
            yticks([0,0.2,0.4,0.6,0.8,1])
            xlim([0,48])
            xticks([0 12 24 36 48])

            hold on
           
          
            

                        
end
hold on
subplot(2,2,2); plot(T, Y(:,28), 'LineWidth', 2, 'Color', 'k')
legend([params{4}(Ym_pert), 'no inhibition'], 'Location', 'SouthEast')
subplot(2,2,4); plot(T, Y(:,28), 'LineWidth', 2, 'Color', 'k')

for v = 1:4
    subplot(2,2,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 12)
end

figure(15)
hold on
for j = 1:length(W_pert)
            dp1 = params;
            dp1{1}(1,W_pert(j)) = dp1{1}(1,W_pert(j))*(1 - 15*1e-2); % perturb each parameter by a small amount

            options = [];
            [dt, dy] = ode23s(@networkODE,tspan,y0,options,dp1);
            dy = real(dy);
            
            hold on
            subplot(2,2,1)
            plot(dt, dy(:,Y_spec(1)), linest(j), 'LineWidth', 2), xlabel('Time (hour)'),  ylabel(append('Normalized ',params{4}(Y_spec(1))))
            set(gca,'FontName','Arial','FontSize',16);
            ylim([0,1.05])
            yticks([0,0.2,0.4,0.6,0.8,1])
            xlim([0,48])
            xticks([0 12 24 36 48])

            dp2 = params;
            dp2{1}(1,W_pert(j)) = dp2{1}(1,W_pert(j))*(1 - 50*1e-2); % perturb each parameter by a small amount

            options = [];
            [dt, dy] = ode23s(@networkODE,tspan,y0,options,dp2);
            dy = real(dy);
            
            hold on
            subplot(2,2,3)
            plot(dt, dy(:,Y_spec(1)), linest(j), 'LineWidth', 2), xlabel('Time (hour)'),  ylabel(append('Normalized ',params{4}(Y_spec(1))))
            set(gca,'FontName','Arial','FontSize',16);
            ylim([0,1.05])
            yticks([0,0.2,0.4,0.6,0.8,1])
            xlim([0,48])
            xticks([0 12 24 36 48])
            hold on
            
                        
end
hold on


subplot(2,2,1); plot(T, Y(:,28), 'LineWidth', 2, 'Color', 'k')
legend([reac_names(W_pert), 'no inhibition'], 'Location', 'SouthEast')
subplot(2,2,3); plot(T, Y(:,28), 'LineWidth', 2, 'Color', 'k')


filename = 'Fig15';
widthInches = 10;
ScriptForExportingImagesForAJP
end



    
