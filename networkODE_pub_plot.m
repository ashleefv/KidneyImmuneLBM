function [Time, Y_pred] = networkODE_pub_plot(tspan, y0, params, tau_index, k_index, n_index, W_index)

GLU = params{1}(1,1);
LPS = params{1}(1,2);

Time = [0:1:48]; % time span

speciesNames = params{4};
%%% reac_names = {'=> GLU','=> LPS','LPS => TLR4','GLU => AGE', 'AGE => RAGE','RAGE => NADPH','TLR4 & ROS => NF\kappaB', 'TLR4 => PI3K','NADPH => ROS', 'PI3K => AKT', 'PI3K => ROS', 'NF\kappaB_e_c => TNF-\alpha','AKT => NF\kappaB','NF\kappaB => IL-6','NF\kappaB => TNF-\alpha','NF\kappaB => VEGF-A_m_R_N_A', 'VEGF-A_m_R_N_A => VEGF-A','NF\kappaB => IL-1\beta','VEGF-A => VEGFR1','VEGF-A => VEGFR2','AGE => RAGE_e_c','RAGE_e_c => NADPH_e_c','VEGFR2 => PI3K_e_c','VEGFR1 => PI3K_e_c', 'NADPH_e_c => ROS_e_c', 'PI3K_e_c => AKT_e_c', 'AKT_e_c => eNOS' , 'VEGFR1 => PLC-\gamma' ,'PLC-\gamma => NF\kappaB_e_c' , 'ROS_e_c => NF\kappaB_e_c','NF\kappaB_e_c => IL-6', 'NF\kappaB_e_c => IL-1\beta', 'eNOS  => NO', 'eNOS => ROS_e_c' ,'ROS_e_c & NO => ONOO','!NO => Ca','PLC-\gamma => Ca', 'Ca => pJunction','pJunction => Gap Width', 'Ca => NO'};
%reac_names = {'$\Rightarrow$ GLU','$\Rightarrow$ LPS','LPS $\Rightarrow$ TLR4','GLU $\Rightarrow$ AGE', 'AGE $\Rightarrow$ RAGE','RAGE $\Rightarrow$ NADPH','TLR4 \& ROS $\Rightarrow$ NF$\kappa$B', 'TLR4 $\Rightarrow$ PI3K','NADPH $\Rightarrow$ ROS', 'PI3K $\Rightarrow$ AKT', 'PI3K $\Rightarrow$ ROS', 'NF$\kappa$B$_{ec}$ $\Rightarrow$ TNF-$\alpha$','AKT $\Rightarrow$ NF$\kappa$B','NF$\kappa$B $\Rightarrow$ IL-6','NF$\kappa$B $\Rightarrow$ TNF-$\alpha$','NF$\kappa$B $\Rightarrow$ VEGF-A$_{mRNA}$', 'VEGF-A$_{mRNA}$ $\Rightarrow$ VEGF-A','NF$\kappa$B $\Rightarrow$ IL-1$\beta$','VEGF-A $\Rightarrow$ VEGFR1','VEGF-A $\Rightarrow$ VEGFR2','AGE $\Rightarrow$ RAGE$_{ec}$','RAGE$_{ec}$ $\Rightarrow$ NADPH$_{ec}$','VEGFR2 $\Rightarrow$ PI3K$_{ec}$','VEGFR1 $\Rightarrow$ PI3K$_{ec}$', 'NADPH$_{ec}$ $\Rightarrow$ ROS$_{ec}$', 'PI3K$_{ec}$ $\Rightarrow$ AKT$_{ec}$', 'AKT$_{ec}$ $\Rightarrow$ eNOS' , 'VEGFR1 $\Rightarrow$ PLC-$\gamma$' ,'PLC-$\gamma$ $\Rightarrow$ NF$\kappa$B$_{ec}$' , 'ROS$_{ec}$ $\Rightarrow$ NF$\kappa$B$_{ec}$','NF$\kappa$B$_{ec}$ $\Rightarrow$ IL-6', 'NF$\kappa$B$_{ec}$ $\Rightarrow$ IL-1$\beta$', 'eNOS  $\Rightarrow$ NO', 'eNOS $\Rightarrow$ ROS$_{ec}$' ,'ROS$_{ec}$ \& NO $\Rightarrow$ ONOO','!NO $\Rightarrow$ Ca','PLC-$\gamma$ $\Rightarrow$ Ca', 'Ca $\Rightarrow$ pJunction','pJunction $\Rightarrow$ Gap Width', 'Ca $\Rightarrow$ NO'};
reac_names = {'\Rightarrow GLU','\Rightarrow LPS','LPS \Rightarrow TLR4','GLU \Rightarrow AGE', 'AGE \Rightarrow RAGE','RAGE \Rightarrow NADPH', 'TLR4 & ROS \Rightarrow NF\kappaB', 'TLR4 \Rightarrow PI3K','NADPH \Rightarrow ROS', 'PI3K \Rightarrow AKT', 'PI3K \Rightarrow ROS', 'NF\kappaB_e_c \Rightarrow TNF-\alpha','AKT \Rightarrow NF\kappaB','NF\kappaB \Rightarrow IL-6','NF\kappaB \Rightarrow TNF-\alpha','NF\kappaB \Rightarrow VEGF-A_m_R_N_A', 'VEGF-A_m_R_N_A \Rightarrow VEGF-A','NF\kappaB \Rightarrow IL-1\beta','VEGF-A \Rightarrow VEGFR1','VEGF-A \Rightarrow VEGFR2','AGE \Rightarrow RAGE_e_c','RAGE_e_c \Rightarrow NADPH_e_c','VEGFR2 \Rightarrow PI3K_e_c','VEGFR1 \Rightarrow PI3K_e_c', 'NADPH_e_c \Rightarrow ROS_e_c', 'PI3K_e_c \Rightarrow AKT_e_c', 'AKT_e_c \Rightarrow eNOS' , 'VEGFR1 \Rightarrow PLC-\gamma' ,'PLC-\gamma \Rightarrow NF\kappaB_e_c' , 'ROS_e_c \Rightarrow NF\kappaB_e_c','NF\kappaB_e_c \Rightarrow IL-6', 'NF\kappaB_e_c \Rightarrow IL-1\beta', 'eNOS  \Rightarrow NO', 'eNOS \Rightarrow ROS_e_c' ,'ROS_e_c & NO \Rightarrow ONOO','!NO \Rightarrow Ca','PLC-\gamma \Rightarrow Ca', 'Ca \Rightarrow pJunction','pJunction \Rightarrow Gap Width', 'Ca \Rightarrow NO'};

size_tau = size(tau_index,2);
size_n = size(n_index,2);
size_k = size(k_index,2); % EC50
size_W = size(W_index, 2);
%% Fitted plots

[Time, Y_pred] = networkODE_run(tspan, y0, params, tau_index, k_index, n_index, W_index, 1);

%% Validation plots

[Time, Y_pred] = networkODE_run(tspan, y0, params, tau_index, k_index, n_index, W_index, 2);

%% Regulatory plots

[Time, Y_pred] = networkODE_run(tspan, y0, params, tau_index, k_index, n_index, W_index, 3);

%% Prediction plots (Fig. 13, 14)

linest = ["-", "--", "-.", ":", "--."];
R1  = [0,0.25,0.5,0.75,1];

if GLU>0 && LPS==0

options = [];
    params{1}(1,1) = 0; params{1}(1,2)=0;
    Y0new(:) = y0; Y0new(:) = 0;
    params{2}(:) = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ]; 

    [To, Yo] = ode23s(@networkODE,tspan,Y0new,options,params);

    figure(13)
    xlim([0,49])
    
    subplot(2,2,1)
    plot(To, Yo(:,23),  linest(1), 'LineWidth', 2)
    hold on
    ylabel("IL-6 activity"); xlabel('Time(hour)')
    set(gca,'FontName','Arial','FontSize',18);
    xticks([0,12,24,36,48])
    ylim([0,1.05])
    grid off;
    hold on

    subplot(2,2,2)
    plot(To, Yo(:,13),  linest(1), 'LineWidth', 2)
    ylabel("ROS_e_c activity"); xlabel('Time(hour)')
    set(gca,'FontName','Arial','FontSize',18);
    xticks([0,12,24,36,48])
    ylim([0,1.05])
    grid off;
    hold on
    
    subplot(2,2,3)
    plot(To, Yo(:,26),  linest(1), 'LineWidth', 2)
    hold on
    ylabel("PLC-\gamma activity"); xlabel('Time(hour)')
    set(gca,'FontName','Arial','FontSize',18);
    xticks([0,12,24,36,48])
    ylim([0,1.05])
    grid off;
    hold on

    subplot(2,2,4)
    plot(To, Yo(:,28),  linest(1), 'LineWidth', 2);  % edited 30 -> 28
    ylabel("pJunction activity"); xlabel('Time(hour)') % edited
    set(gca,'FontName','Arial','FontSize',18);
    xticks([0,12,24,36,48]); ylim([0,1.05]);
    legendI{1} = ['W_G_L_U = ' num2str(R1(1))];

for i = [2:length(R1)]
     
    % tau values for GLU only (reset to nominal values)
    params{2}(:) = [1	5.06	1	1	1	3.39	1	1	5.04	1	1	9.97	5.89	5.04	1	6.60	5.75	1	5.07	4.36	1	1.91	4.66	9.90	4.90	1	3.14	1	1	1];

    
    
    params{1}(1,1) = R1(i);
    params{1}(1,2) = 0;
  
    
    options = [];
    [To, Yo] = ode23s(@networkODE,tspan,y0,options,params);

    figure(13)
    xlim([0,49])
    
    hold on
    subplot(2,2,1)
    plot(To, Yo(:,23),  linest(i), 'LineWidth', 2)
    hold on
    ylabel("IL-6 activity"); xlabel('Time (hour)')
    set(gca,'FontName','Arial','FontSize',18);
    xticks([0,12,24,36,48]); xlim([0,49]); ylim([0,1.05])

    grid off;
    hold on

    subplot(2,2,2)
    hold on
    plot(To, Yo(:,13),  linest(i), 'LineWidth', 2)
    ylabel("ROS_e_c activity"); xlabel('Time (hour)')
    set(gca,'FontName','Arial','FontSize',18);
    xticks([0,12,24,36,48]); xlim([0,49]); ylim([0,1.05])

    grid off;
    hold on
    
    subplot(2,2,3)
    hold on
    plot(To, Yo(:,26),  linest(i), 'LineWidth', 2)
    hold on
    ylabel("PLC-\gamma activity"); xlabel('Time (hour)')
    set(gca,'FontName','Arial','FontSize',18);
    xticks([0,12,24,36,48]); xlim([0,49]); ylim([0,1.05])

    grid off;
    hold on

    subplot(2,2,4)
    hold on
    plot(To, Yo(:,28),  linest(i), 'LineWidth', 2)
    ylabel("pJunction activity"); xlabel('Time (hour)')
    set(gca,'FontName','Arial','FontSize',18);
    xticks([0,12,24,36,48]); xlim([0,49]); ylim([0,1.05])

    grid off;
    hold on
    
    legendI{i} = ['W_G_L_U = ' num2str(R1(i))];
    legend(legendI)
    hold on
end

for i = [1:length(R1)]

    params{1}(1,1) = 1;
    params{1}(1,2) = R1(i);
    
    options = [];
    [To, Yo] = ode23s(@networkODE,tspan,y0,options,params);

    figure(14)

    subplot(2,2,1)
    plot(To, Yo(:,24),  linest(i), 'LineWidth', 2)
    hold on
    ylabel("TNF-\alpha activity"); xlabel('Time (hour)')
    set(gca,'FontName','Arial','FontSize',18);
    xticks([0,12,24,36,48]);     xlim([0,49])


    grid off;
    hold on

    subplot(2,2,2)
    plot(To, Yo(:,25),  linest(i), 'LineWidth', 2)
    ylabel("IL-1\beta activity"); xlabel('Time (hour)')
    set(gca,'FontName','Arial','FontSize',18);
    xticks([0,12,24,36,48]);    xlim([0,49])

    grid off;
    hold on
    
    subplot(2,2,3)
    plot(To, Yo(:,26),  linest(i), 'LineWidth', 2)
    hold on
    ylabel("PLC-\gamma activity"); xlabel('Time (hour)')
    set(gca,'FontName','Arial','FontSize',18);
    xticks([0,12,24,36,48]);    xlim([0,49])


    grid off;
    hold on

    subplot(2,2,4)
    plot(To, Yo(:,20),  linest(i), 'LineWidth', 2)
    ylabel("NO activity"); xlabel('Time (hour)')
    set(gca,'FontName','Arial','FontSize',18);
    xticks([0,12,24,36,48]);     xlim([0,49])


    grid off;
    hold on
    legendInfo{i} = ['W_{GLU}, W_{LPS}=(1,' num2str(R1(i)) ')'];

    hold on

end
hold on
legend(legendInfo)

end
%% Prediction plots (Fig. 15)
if GLU>0 && LPS>0
    sens_change = -15;

    [s_FD_Ym, s_FD_W, s_FD_k, W_pert, Ym_pert] = post_sens(params, y0, tspan, sens_change);
end
%% UQLab plots
% 8 plots

UQLab_plt = 'ON';
if strcmp(UQLab_plt, 'ON')
    disp('(1) Install UQLab.')
    disp('(2) addpath("~/UQLab_Rel2.0.0/UQLab_Rel2.0.0/core")')
    disp('(3) Run UQLab_Sobolplots_networkmodel() in UQLab_scripts folder.')
end

%% Histograms (FIG )

if GLU>0 && LPS==0
    load data/pub_plot_GLU.mat
    accep_mean = mean(p_posterior(:,[1:end]));

elseif GLU==0 && LPS>0
    load data/pub_plot_LPS.mat
    accep_mean = mean(p_posterior(:,[1:end]));

elseif GLU>0 && LPS>0
    load data/pub_plot_both.mat
    accep_mean = mean(p_posterior(:,[1:end]));
else
    disp('"insufficient initialization", choose a treatment condition.');
end

figure(54) % fig S4

hold on
for i=1:size_tau
    subplot(3,5,i)
    grid on
    hold on
    % hist(population(:,i)) % outdated histogram function plots with navy blue - color
    histogram(population(:,i), 'FaceColor', '#020259');
    hold on
    xline(accep_mean(:,i), 'LineWidth', 2, 'color', 'r')
    hold on
    title([num2str(tau_index(i)) + ": " + params{4}(tau_index(i))])
    hold on
end
hold on
legend('\tau', 'mean', 'FontSize', 10)

if GLU>0 && LPS==0

    figure(531) % Fig S3 (1)
    hold on
    for i=1:size_W
        subplot(3,5,i)
        grid on
        hold on
        histogram(population(:,i+size_tau), 'FaceColor', '#00FFFF')
        hold on
        xline(accep_mean(:,i+size_tau), 'LineWidth', 2, 'color', 'r')
        hold on
        title([num2str(W_index(i)) + ": " + reac_names(W_index(i))])
    end
    hold on
    legend('W', 'mean', 'FontSize', 10)


    figure(532) % Fig S3 (2)
    hold on
    for i=1:size_n
        subplot(4,8,i)
        grid on
        hold on
        histogram(population(:,i+size_tau+size_W), 'FaceColor', '#EDB120')
        hold on
        xline(accep_mean(i+size_tau+size_W), 'LineWidth', 2, 'color', 'r')
        hold on
        title([num2str(n_index(i)) + ": " + reac_names(n_index(i))])
    end
    hold on
    legend('n', 'mean', 'FontSize', 10)

    figure(533) % fig S3 (3)  
    hold on
    for i=1:size_k
        subplot(4,6,i)
        grid on
        hold on
        histogram(population(:,i+size_tau+size_W+size_n), 'FaceColor', '#77AC30')
        hold on
        xline(accep_mean(i+size_tau+size_W+size_n), 'LineWidth', 2, 'color', 'r')
        hold on
        title([num2str(k_index(i)) + ": " + reac_names(k_index(i))])
    end
    hold on
    legend('EC_{50}', 'mean', 'FontSize', 10)

end

% posterior prediction at 12 hours

figure(55) % Fig S5
vars = [23, 24, 25, 6, 13, 22, 20, 12, 27, 28]; % edited 30 -> 28
hold on
for i=1:10
    subplot(2,5,i)
    grid on
    hold on
    histogram(abs(Yp(:,12,vars(i))), 'FaceColor', '#7E2F8E')
    hold on
    title([num2str(vars(i)) + ": " + params{4}(vars(i))])
end
hold on
legend('prediction at 12 hours', 'FontSize', 10)


end
