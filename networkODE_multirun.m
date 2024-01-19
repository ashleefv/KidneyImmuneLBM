function [p_posterior,population, Yp, credible] = networkODE_multirun(tspan, y0, params, p_fitted, error_fitted, global_p_best, tau_index, k_index, n_index, W_index, mode)

load invitro_data.mat

GLU = params{1}(1,1);
LPS = params{1}(1,2);

Time = [0:1:48]; % time span

speciesNames = params{4};
% reac_names = {'=> GLU','=> LPS','LPS => TLR4','GLU => AGE', 'AGE => RAGE','RAGE => NADPH','TLR4 & ROS => NF\kappaB', 'TLR4 => PI3K','NADPH => ROS', 'PI3K => AKT', 'PI3K => ROS', 'NF\kappaB_e_c => TNF-\alpha','AKT => NF\kappaB','NF\kappaB => IL-6','NF\kappaB => TNF-\alpha','NF\kappaB => VEGF-A_m_R_N_A', 'VEGF-A_m_R_N_A => VEGF-A','NF\kappaB => IL-1\beta','VEGF-A => VEGFR1','VEGF-A => VEGFR2','AGE => RAGE_e_c','RAGE_e_c => NADPH_e_c','VEGFR2 => PI3K_e_c','VEGFR1 => PI3K_e_c', 'NADPH_e_c => ROS_e_c', 'PI3K_e_c => AKT_e_c', 'AKT_e_c => eNOS' , 'VEGFR1 => PLC-\gamma' ,'PLC-\gamma => NF\kappaB_e_c' , 'ROS_e_c => NF\kappaB_e_c','NF\kappaB_e_c => IL-6', 'NF\kappaB_e_c => IL-1\beta', 'eNOS  => NO', 'eNOS => ROS_e_c' ,'ROS_e_c & NO => ONOO','!NO => Ca','PLC-\gamma => Ca', 'Ca => pJunction','pJunction => Gap Width', 'Ca => NO'};
reac_names = {'\Rightarrow GLU','\Rightarrow LPS','LPS \Rightarrow TLR4','GLU \Rightarrow AGE', 'AGE \Rightarrow RAGE','RAGE \Rightarrow NADPH', 'TLR4 & ROS \Rightarrow NF\kappaB', 'TLR4 \Rightarrow PI3K','NADPH \Rightarrow ROS', 'PI3K \Rightarrow AKT', 'PI3K \Rightarrow ROS', 'NF\kappaB_e_c \Rightarrow TNF-\alpha','AKT \Rightarrow NF\kappaB','NF\kappaB \Rightarrow IL-6','NF\kappaB \Rightarrow TNF-\alpha','NF\kappaB \Rightarrow VEGF-A_m_R_N_A', 'VEGF-A_m_R_N_A \Rightarrow VEGF-A','NF\kappaB \Rightarrow IL-1\beta','VEGF-A \Rightarrow VEGFR1','VEGF-A \Rightarrow VEGFR2','AGE \Rightarrow RAGE_e_c','RAGE_e_c \Rightarrow NADPH_e_c','VEGFR2 \Rightarrow PI3K_e_c','VEGFR1 \Rightarrow PI3K_e_c', 'NADPH_e_c \Rightarrow ROS_e_c', 'PI3K_e_c \Rightarrow AKT_e_c', 'AKT_e_c \Rightarrow eNOS' , 'VEGFR1 \Rightarrow PLC-\gamma' ,'PLC-\gamma \Rightarrow NF\kappaB_e_c' , 'ROS_e_c \Rightarrow NF\kappaB_e_c','NF\kappaB_e_c \Rightarrow IL-6', 'NF\kappaB_e_c \Rightarrow IL-1\beta', 'eNOS  \Rightarrow NO', 'eNOS \Rightarrow ROS_e_c' ,'ROS_e_c & NO \Rightarrow ONOO','!NO \Rightarrow Ca','PLC-\gamma \Rightarrow Ca', 'Ca \Rightarrow pJunction','pJunction \Rightarrow Gap Width', 'Ca \Rightarrow NO'};

size_tau = size(tau_index,2);
size_n = size(n_index,2);
size_k = size(k_index,2); % EC50
size_W = size(W_index, 2);


if GLU>0 && LPS==0                                   % condition (GLU = 1, LPS = 0)
    
    t_data = train_data([2,6,10,14,18,22],:);        % Training data for given time points
    e_data = Terror([2,6,10,14,18,22],:);            % Error bars for training
    v_data = v_data([2,6,10,14,18],:);               % Validation data for given time points
    ev_data = Verror([2,6,10,14,18],:);              % Error bars for validation
    timedata = time_data([2,6,10,14,18,22],:);       % Time points for training set
    timev = vtime([2,6,10,14,18],:);                 % Time points for validation set

end 

if GLU==0 && LPS>0                                   % condition (LPS = 1, GLU = 0)
   
    t_data = train_data([3,7,11,15,19,23],:);        % Training data for given time points
    e_data = Terror([3,7,11,15,19,23],:);            % Error bars for training
    v_data = v_data([3,7,11,15,19],:);               % Validation data for given time points
    ev_data = Verror([3,7,11,15,19],:);              % Error bars for validation
    timedata = time_data([3,7,11,15,19,23],:);       % Time points for training set
    timev = vtime([3,7,11,15,19],:);                 % Time points for validation set

end

if GLU>0 && LPS>0                                    % condition (GLU = 1, LPS = 1)
    t_data = train_data([4,8,12,16,20,24],:);        % Training data for given time points
    e_data = Terror([4,8,12,16,20,24],:);            % Error bars for training
    v_data = v_data([4,8,12,16,20],:);               % Validation data for given time points
    ev_data = Verror([4,8,12,16,20],:);              % Error bars for validation
    timedata = time_data([4,8,12,16,20,24],:);       % Time points for training set
    timev = vtime([4,8,12,16,20],:);                 % Time points for validation set

end


options=[];
[T_best, Y_best] = ode23s(@networkODE,tspan,y0,options,params);

%%
if GLU>0 && LPS==0
    GB = params;
    GB{2}(tau_index) = global_p_best(1:size_tau);
    GB{1}(1,W_index) = global_p_best(size_tau+1:size_tau+size_W);
    GB{1}(2,n_index) = global_p_best(size_tau+size_W+1:size_tau+size_W+size_n);
    GB{1}(3,k_index) = global_p_best(size_tau+size_W+size_n+1:size_tau+size_W+size_n+size_k);
end
if GLU==0 && LPS>0
    GB = params;
    GB{2}(tau_index) = global_p_best(1:size_tau);
end

if GLU>0 && LPS>0
    GB = params;
    GB{2}(tau_index) = global_p_best(1:size_tau);
    
end

% sort parameters by minimum error
col = [0 0.4470 0.7410];
Params_error = [error_fitted p_fitted]; % coef_fitted is the parameter output from all the LHS fmincon runs
error_column = 1; % first column stores error_fitted
[sorted_Params_error, index] = sortrows(Params_error,error_column);
param_posterior = sorted_Params_error(:,[2:end]);

options = [];
[T_gb, Y_gb] = ode23s(@networkODE,tspan,y0,options,GB); % prediction at global-best estimates


% Sum of Squared Error (SSE) plot
Params_error = [error_fitted p_fitted]; % coef_fitted is the parameter output from all the LHS fmincon runs
error_column = 1; % first column stores error_fitted
[sorted_Params_error, index] = sortrows(Params_error,error_column);
accep_id = find(sorted_Params_error(:,1) < 1.2*min(sorted_Params_error(:,1)));

figure(20)
%-- generate figure to show the error of all the parameter values from LHS
for i = 1:size(p_fitted, 2)
    plot(1:length(index),sorted_Params_error(:,error_column),'o')
    hold on
    xlabel('number of runs')
    ylabel('sum of squared error (sse)')
end

% Monte Carlo Method
disp('Starting Monte Carlo Simulations ...')

p_posterior = sorted_Params_error(accep_id,[2:end]); % acceptable parameters within 1.2*min(SSE)
Ns = 5000; % number of samples
population = zeros(1,length(global_p_best));         % initialize population
%population = zeros(1,size(p_fitted,2));
s = RandStream('mlfg6331_64');
tic
for j = 1:Ns
    for m = 1:length(global_p_best)
        population(j,m) = randsample(s, p_posterior(:,m), 1); % posterior population of acceptable parameters
    end
end
sprintf('population size: %d x %d',size(population))

for j = 1:Ns
    if GLU>0 && LPS==0
        dp = params;
        dp{2}(tau_index) = population(j,(1:size_tau));
        dp{1}(1,W_index) = population(j,(size_tau+1:size_tau+size_W));
        dp{1}(2,n_index) = population(j,(size_tau+size_W+1:size_tau+size_W+size_n));
        dp{1}(3,k_index) = population(j,(size_tau+size_W+size_n+1:end)); 

        options = [];
        [Tout, Yout] = ode23s(@networkODE,tspan,y0,options,dp);
        Yp(j,:,:) = Yout;
        dp = params;
    end
    if GLU==0 && LPS>0
        dp = params;
        dp{2}(tau_index) = population(j,(1:size_tau));
        
        options = [];
        [Tout, Yout] = ode23s(@networkODE,tspan,y0,options,dp);
        Yp(j,:,:) = Yout;
        dp = params;
    end
    if GLU>0 && LPS>0
        dp = params;
        dp{2}(tau_index) = population(j,(1:size_tau));
        
        options = [];
        [Tout, Yout] = ode23s(@networkODE,tspan,y0,options,dp);
        Yp(j,:,:) = Yout;
        dp = params;
    end
    
    
   
end
toc
%%
% Time = [0:1:48];
% vars = [23, 24, 25, 6, 13, 22, 20, 12, 27, 28];

% for v = 1:length(vars)
  %  for i = 1:length(Time)
    
   %     mean_y(i,v) = mean(abs(Yp(:,i,vars(v))));
   %     sd_y(i,v) = std(abs(Yp(:,i,vars(v))));
        

   % end
   % ts = tinv([0.025  0.975],length(Time)-1);
   % CI95 = abs(mean_y(:,v) + ts.*(sd_y(:,v)/sqrt(length(Time))));
   % lower_CI(:,v) = CI95(:,1);
   % upper_CI(:,v) = CI95(:,2);
% end

         
% save('workspaces/both_best_fitted_MC')
%%
vars = [23, 24, 25, 6, 13, 22, 20, 12, 27, 28];

blue = 	[0 0.4470 0.7410];
for v = 1:8
    for t = 1:length(Time)
    
       credible(t,:) = quantile(abs(Yp(:,t,vars(v))), [0.025 0.975]);
    end
    
    figure(21)
    grid on
    hold on
    subplot(2,4,v)
    hold on
    % --- Prediction
    % plot(Time, mean_y(:,v), 'color', 'b', 'LineWidth', 0.8)            % mean of prediction posterior
    % hold on
    plot(T_gb, Y_gb(:,vars(v)), '--', 'color', 'b', 'LineWidth', 0.8)   % global best estimates
    hold on
    plot(T_best, Y_best(:,v), 'color', 'k', 'LineWidth', 1)             % mean of acceptable estimates
    hold on
    [ph,msg] = jbfill(Time,credible(:,1)', credible(:,2)',blue,blue,1,0.2);
    hold on    

end

%% Histograms

accep_mean = mean(p_posterior(:,[1:end]));

figure(22)

hold on
for i=1:size_tau
subplot(3,5,i)
grid on
hold on
histogram(population(:,i), 'FaceColor', '#020259')
hold on
xline(accep_mean(:,i), 'LineWidth', 2, 'color', 'r')
hold on
title([num2str(tau_index(i)) + ": " + params{4}(tau_index(i))])
hold on
end
hold on
legend('\tau', 'mean', 'FontSize', 10)

if GLU>0 && LPS==0

figure(23)
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


figure(24)
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


figure(25)
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

figure(26)  
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
%% posterior prediction at 12 hours

% figure(27)
% vars = [23, 24, 25, 6, 13, 22, 20, 12, 27, 28];
% hold on
%for i=1:9
%    subplot(3,3,i)
%    grid on
%    hold on
%    histogram(abs(Yp(:,12,vars(i))), 'FaceColor', '#7E2F8E')
%    hold on
%    title([num2str(vars(i)) + ": " + params{4}(vars(i))])
%end
%hold on
%legend('prediction at 12 hours', 'FontSize', 10)



%%

if mode==1

            figure(28)
            hold on
            vars = [23, 24, 25, 6, 13, 22, 20, 12, 27, 28];

             subplot(2,4,1)
             hold on
             % --- Prediction
             plot(T_best, Y_best(:,23), 'color', 'k',  'LineWidth', 1.2)
             hold on
             % --- Training data
             scatter(timedata(:,1), t_data(:,1), 50, 'filled')
             yerr = abs(e_data(:,1)-t_data(:,1));
             errorbar(timedata(:,1), t_data(:,1), yerr, 'o');
             hold on            
             xlabel('Time (hour)');
             ylabel(append(speciesNames([23]), ' activity'));
             grid off
             ylim([0,1.05])
             xlim([0,48])
             xticks([0 12 24 36 48])
             set(gca,'FontName','Arial','FontSize',18);
             t = title('Tessaro et al. (2020)'); t.FontSize = 14;
             
                       
             subplot(2,4,2)
             hold on
             % --- Prediction
             plot(T_best, Y_best(:,24) ,'color', 'k', 'LineWidth', 1.2)
             hold on
             % --- Training data
             scatter(timedata(:,3), t_data(:,3),  50, 'filled')
             yerr = abs(e_data(:,3)-t_data(:,3));
             errorbar(timedata(:,3), t_data(:,3),  yerr, 'o');
             hold on
             grid off
             ylim([0,1.05])
             xlim([0,48])
             xticks([0 12 24 36 48])
             xlabel('Time (hour)');
             ylabel(append(speciesNames([24]), ' activity'));
             set(gca,'FontName','Arial','FontSize',18);
             t = title('Tessaro et al. (2020)'); t.FontSize = 14;
             
          
     
             subplot(2,4,3)
             hold on
             % --- Prediction
             plot(T_best, Y_best(:,25),  'color', 'k', 'LineWidth', 1.2)
             hold on
             % --- Training data
             scatter(timedata(:,2), t_data(:,2),  50, 'filled')
             yerr = abs(e_data(:,2)-t_data(:,2));
             errorbar(timedata(:,2), t_data(:,2), yerr, 'o');
             hold on
             grid off
             ylim([0,1.05])
             xlim([0,48])
             xticks([0 12 24 36 48])
             xlabel('Time (hour)');
             ylabel(append(speciesNames([25]), ' activity')); 
             set(gca,'FontName','Arial','FontSize',18); 
             t = title('Li et al. (2021)'); t.FontSize = 14;
            

         
             subplot(2,4,4)
             hold on
             % --- Prediction
             plot(T_best, Y_best(:,6),  'color', 'k', 'LineWidth', 1.2)
             hold on
             % --- Training data
             scatter(timedata(:,4), t_data(:,4), 50, 'filled')
             yerr = abs(e_data(:,4)-t_data(:,4));
             errorbar(timedata(:,4), t_data(:,4), yerr, 'o');
             hold on
             grid off
             ylim([0,1.05])
             xlim([0,48])
             xticks([0 12 24 36 48])
             xlabel('Time (hour)');
             ylabel(append(speciesNames([6]), ' activity')); 
             set(gca,'FontName','Arial','FontSize',18); 
             t = title('Li et al. (2021)'); t.FontSize = 14;
            

             subplot(2,4,5);
             % --- Prediction
             hold on
             plot(T_best, Y_best(:,13),  'color', 'k', 'LineWidth', 1.2)
             hold on
             % --- Training data
             scatter(timedata(:,5), t_data(:,5),  50, 'filled')
             yerr = abs(e_data(:,5)-t_data(:,5));
             errorbar(timedata(:,5), t_data(:,5), yerr, 'o');
             hold on 
             grid off
             ylim([0,1.05])
             xlim([0,48])
             xticks([0 12 24 36 48])
             xlabel('Time (hour)');
             ylabel(append(speciesNames([13]), ' activity')); 
             set(gca,'FontName','Arial','FontSize',18);
             t = title('de Souza et al. (2007)'); t.FontSize = 14;
             

             subplot(2,4,6);
             % --- Prediction
             hold on
             plot(T_best, Y_best(:,22),   'color', 'k', 'LineWidth', 1.2)
             hold on
             % --- Training data
             scatter(timedata(:,6), t_data(:,6),  50, 'filled')
             yerr = abs(e_data(:,6)-t_data(:,6));
             errorbar(timedata(:,6), t_data(:,6), yerr, 'o');
             hold on
             grid off
             ylim([0,1.05])
             xlim([0,48])
             xticks([0 12 24 36 48])
             xlabel('Time (hour)');
             ylabel(append(speciesNames([22]), ' activity')); 
             set(gca,'FontName','Arial','FontSize',18);
             t = title('Nakagawa et al. (2006)'); t.FontSize = 14;

           

              subplot(2,4,7)
              hold on
              % --- Prediction
              plot(T_best, Y_best(:,20),   'color', 'k', 'LineWidth', 1.2)
              hold on
              % --- Training data
              scatter(timedata(:,7),  t_data(:,7), 50, 'filled')
              yerr = abs(e_data(:,7)-t_data(:,7));
              errorbar(timedata(:,7),  t_data(:,7), yerr, 'o');
              hold on
              grid off
              ylim([0,1.05])
              xlim([0,48])
              xticks([0 12 24 36 48])
              xlabel('Time (hour)');
              ylabel(append(speciesNames([20]), ' activity')); 
              set(gca,'FontName','Arial','FontSize',18);
              t = title('Kasai et al. (2005)'); t.FontSize = 14;
              grid off

              subplot(2,4,8);
              hold on
              % --- Prediction
              plot(T_best, Y_best(:,12),  'color', 'k', 'LineWidth', 1.2)
              hold on
              % --- Training data
              scatter(timedata(:,9),  t_data(:,9), 50, 'filled')
              yerr = abs(e_data(:,9)-t_data(:,9));
              errorbar(timedata(:,9),  t_data(:,9), yerr, 'o');
              hold on
              grid off
              ylim([0,1.05])
              xlim([0,48])
              xticks([0 12 24 36 48])
              xlabel('Time (hour)');
              ylabel(append(speciesNames([12]), ' activity'));
              set(gca,'FontName','Arial','FontSize',18);  
              hold on
              t = title('Yang et al. (2008)'); t.FontSize = 14;

              blue = [0 0.4470 0.7410];
              for v = 1:8
                  for t = 1:length(Time)
                      credible(t,:) = quantile(abs(Yp(:,t,vars(v))), [0.025 0.975]);
                  end
                  figure(28)
                  hold on
                  subplot(2,4,v)
                  hold on
                  [ph,msg] = jbfill(Time,credible(:,1)', credible(:,2)',blue,blue,1,0.2);
                  hold on    

              end
              lgd = legend('model prediction', 'data', 'error bar', '95% credible interval',  'Location', 'SouthEast');
              lgd.FontSize = 12;

elseif mode==2
    purp = [0.4940    0.1840    0.5560]; % color

    

    figure(29)

    if GLU>0 && LPS==0

             vars = [23,24,25,27,12,20];

             subplot(2,3,1)
             % --- Prediction
             plot(T_best, Y_best(:,vars(1)), 'color', 'k', 'LineWidth', 1.2)
             hold on
             %plot(Time, mean_y(:,1), 'color', purp, 'LineWidth', 0.8)
             %hold on
             
             % --- Validation data
             scatter(timev(:,1), v_data(:,1), 30, 'ro', 'MarkerFaceColor', 'r')
             yerr = abs(ev_data(:,1) - v_data(:,1));
             errorbar(timev(:,1), v_data(:,1), yerr, 'ro');
             hold on                      
             xlabel('Time (hour)');
             ylabel(append(speciesNames([23]), ' activity'));
             ylim([0,1.05])
              xlim([0,48])
              xticks([0 12 24 36 48])
             set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Ayala et al. (2019)'); t.FontSize = 14;
             grid off
            
             subplot(2,3,2)
             % --- Prediction
             plot(T_best, Y_best(:,vars(2)), 'color', 'k', 'LineWidth', 1.2)
             hold on

             %plot(Time, mean_y(:,2), 'color', purp, 'LineWidth', 0.8)
             %hold on
             % --- Validation data
             scatter(timev(:,3), v_data(:,3),  30, 'ro', 'MarkerFaceColor', 'r')
             yerr = abs(ev_data(:,3) - v_data(:,3));
             errorbar(timev(:,3), v_data(:,3),  yerr, 'ro');
             hold on         
             xlabel('Time (hour)');
             ylabel(append(speciesNames([24]), ' activity'));
             ylim([0,1.05])
              xlim([0,48])
              xticks([0 12 24 36 48])
             set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Ayala et al. (2019)'); t.FontSize = 14;             
             grid off
          
         
             subplot(2,3,3);
             % --- Prediction
             plot(T_best, Y_best(:,vars(3)), 'color', 'k', 'LineWidth', 1.2)
             hold on
             %plot(Time, mean_y(:,3),'color', purp, 'LineWidth', 0.8)
             %hold on
             % --- Validation data
             scatter(timev(:,2), v_data(:,2),  30, 'ro', 'MarkerFaceColor', 'r')
             yerr = abs(ev_data(:,2)-v_data(:,2));
             errorbar(timev(:,2), v_data(:,2), yerr, 'ro');
             hold on
             xlabel('Time (hour)');
             ylabel(append(speciesNames([25]), ' activity')); 
             ylim([0,1.05])
              xlim([0,48])
              xticks([0 12 24 36 48])
              set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('He et al. (2020)'); t.FontSize = 14;
             grid off

         
             subplot(2,3,4);
             % --- Prediction
             plot(T_best, Y_best(:,vars(9)), 'color', 'k', 'LineWidth', 1.2)
             hold on
             %plot(Time, mean_y(:,9), 'color', purp, 'LineWidth', 0.8)
             %hold on
             % --- Validation data
             scatter(timev(:,4), v_data(:,4), 30, 'ro', 'MarkerFaceColor', 'r')
             yerr = abs(ev_data(:,4)-v_data(:,4));
             errorbar(timev(:,4), v_data(:,4), yerr,'ro');
             hold on
             xlabel('Time (hour)');
             ylabel(append(speciesNames([27]), ' activity')); 
             ylim([0,1.05])
              xlim([0,48])
              xticks([0 12 24 36 48])
              set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Yang et al. (2008)'); t.FontSize = 14;
             grid off


             subplot(2,3,5)
             % --- Prediction
             plot(T_best, Y_best(:,vars(8)), 'color', 'k', 'LineWidth', 1.2)
             hold on
             %plot(Time, mean_y(:,8), 'color', purp, 'LineWidth', 0.8)
             %hold on
              % --- Validation data
              scatter(timev(:,5), v_data(:,5), 30, 'ro', 'MarkerFaceColor', 'r')
              yerr = abs(ev_data(:,5)-v_data(:,5));
              errorbar(timev(:,5), v_data(:,5), yerr, 'ro');
              hold on              
              xlabel('Time (hour)');
              ylabel(append(speciesNames([12]), ' activity'));
              ylim([0,1.05])
              xlim([0,48])
              xticks([0 12 24 36 48])
              set(gca,'FontName','Arial','FontSize',18);  
              hold on
              t = title('Seo et al. (2022)'); t.FontSize = 14;
              hold on
            
              grid off

              subplot(2,3,6)
              % --- Prediction
              plot(T_best, Y_best(:,vars(7)), 'color', 'k', 'LineWidth', 1.2)
              hold on             
              %plot(Time, mean_y(:,7), 'color', purp, 'LineWidth', 0.8)
              %hold on
              % --- Validation data
              scatter(timev(:,7),  v_data(:,7), 30, 'ro', 'MarkerFaceColor', 'r')
              yerr = abs(ev_data(:,7)-v_data(:,7));
              errorbar(timev(:,7),  v_data(:,7), yerr, 'ro');
              hold on         
              xlabel('Time (hour)');
              ylabel(append(speciesNames([20]), ' activity')); 
              ylim([0,1.05])
              xlim([0,48])
              xticks([0 12 24 36 48])
              set(gca,'FontName','Arial','FontSize',18);
              hold on
              t = title('Ramanathan et al. (2002)'); t.FontSize = 14;
              
              grid off

            
            for v = 1:length(vars)
                for t = 1:length(Time)
                    credible(t,:) = quantile(abs(Yp(:,t,vars(v))), [0.025 0.975]);
                end
                figure(29)
                grid off
                hold on
                subplot(2,3,v)
                hold on
                [ph,msg] = jbfill(Time,credible(:,1)', credible(:,2)',purp, purp,1,0.2);
                hold on
            end
            lgd = legend('model prediction', 'data', 'error bar', '95% credible interval',  'Location', 'SouthEast');
            lgd.FontSize = 14;
    end
    if GLU==0 && LPS>0

        vars = [23,24,25,27,28,20];

        subplot(2,3,1)
             % --- Prediction
             plot(T_best, Y_best(:,vars(1)), 'color', 'k', 'LineWidth', 1.2)
             hold on

             %plot(Time, mean_y(:,1), 'color', purp, 'LineWidth', 0.8)
             %hold on
             % --- Validation data
             scatter(timev(:,1), v_data(:,1), 30, 'ro', 'MarkerFaceColor', 'r')
             yerr = abs(ev_data(:,1) - v_data(:,1));
             errorbar(timev(:,1), v_data(:,1), yerr, 'ro');
             hold on                      
             xlabel('Time (hour)');
             ylabel(append(speciesNames([23]), ' activity'));
             ylim([0,1.05])
              xlim([0,48])
              xticks([0 12 24 36 48])
             set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Ayala et al. (2019)'); t.FontSize = 14;
             grid off
            
             subplot(2,3,2)
             % --- Prediction
             plot(T_best, Y_best(:,vars(2)),'color', 'k', 'LineWidth', 1.2)
             hold on
             %plot(Time, mean_y(:,2),'color', purp, 'LineWidth', 0.8)
             %hold on
             % --- Validation data
             scatter(timev(:,3), v_data(:,3),  30, 'ro', 'MarkerFaceColor', 'r')
             yerr = abs(ev_data(:,3) - v_data(:,3));
             errorbar(timev(:,3), v_data(:,3),  yerr, 'ro');
             hold on             
            
             xlabel('Time (hour)');
             ylabel(append(speciesNames([24]), ' activity'));
             ylim([0,1.05])
              xlim([0,48])
              xticks([0 12 24 36 48])
              set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Ayala et al. (2019)'); t.FontSize = 14;             
             grid off
          
         
             subplot(2,3,3);
             % --- Prediction
             plot(T_best, Y_best(:,vars(3)), 'color', 'k', 'LineWidth', 1.2)
             hold on

             %plot(Time, mean_y(:,3),'color', purp, 'LineWidth', 0.8)
             %hold on
             % --- Validation data
             scatter(timev(:,2), v_data(:,2),  30, 'ro', 'MarkerFaceColor', 'r')
             yerr = abs(ev_data(:,2)-v_data(:,2));
             errorbar(timev(:,2), v_data(:,2), yerr, 'ro');
             hold on             
             
             xlabel('Time (hour)');
             ylabel(append(speciesNames([25]), ' activity')); 
             ylim([0,1.05])
              xlim([0,48])
              xticks([0 12 24 36 48])
              set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('He et al. (2020)'); t.FontSize = 14;
             grid off

         
             subplot(2,3,4);
             % --- Prediction
             plot(T_best, Y_best(:,vars(4)),'color', 'k', 'LineWidth', 1.2)
             hold on
             %plot(Time, mean_y(:,9),'color', purp, 'LineWidth', 0.8)
             %hold on
             % --- Validation data
             scatter(timev(:,4), v_data(:,4), 30, 'ro', 'MarkerFaceColor', 'r')
             yerr = abs(ev_data(:,4)-v_data(:,4));
             errorbar(timev(:,4), v_data(:,4), yerr,'ro');
             hold on
             xlabel('Time (hour)');
             ylabel(append(speciesNames([27]), ' activity')); 
             ylim([0,1.05])
              xlim([0,48])
              xticks([0 12 24 36 48])
              set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Yang et al. (2008)'); t.FontSize = 14;
             grid off


             subplot(2,3,5)
             % --- Prediction
              plot(T_best, Y_best(:,vars(5)), 'color', 'k', 'LineWidth', 1.2)
              hold on
              %plot(Time, mean_y(:,10), 'color', purp, 'LineWidth', 0.8)
              %hold on   
              % --- Validation data
              scatter(timev(:,8), v_data(:,8), 30, 'ro', 'MarkerFaceColor', 'r')
              yerr = abs(ev_data(:,8)-v_data(:,8));
              errorbar(timev(:,8), v_data(:,8), yerr, 'ro');
              hold on              
              xlabel('Time (hour)');
              ylabel('VE-Cadherin activity');
              ylim([0,1.05])
              xlim([0,48])
              xticks([0 12 24 36 48])
              set(gca,'FontName','Arial','FontSize',18);  
              hold on
              t = title('Seo et al. (2022)'); t.FontSize = 14;
              hold on
            
              grid off

              subplot(2,3,6)
              % --- Prediction
              plot(T_best, Y_best(:,vars(6)), 'color', 'k', 'LineWidth', 1.2)
              hold on             
              %plot(Time, mean_y(:,7), 'color', purp, 'LineWidth', 0.8)
              %hold on
              % --- Validation data
              scatter(timev(:,7),  v_data(:,7), 30, 'ro', 'MarkerFaceColor', 'r')
              yerr = abs(ev_data(:,7)-v_data(:,7));
              errorbar(timev(:,7),  v_data(:,7), yerr, 'ro');
              hold on                            
              xlabel('Time (hour)');
              ylabel(append(speciesNames([20]), ' activity')); 
              ylim([0,1.05])
              xlim([0,48])
              xticks([0 12 24 36 48])
              set(gca,'FontName','Arial','FontSize',18);
              hold on
              t = title('Ramanathan et al. (2002)'); t.FontSize = 14;
              grid off

              
              for v = 1:length(vars)
                for t = 1:length(Time)
                    credible(t,:) = quantile(abs(Yp(:,t,vars(v))), [0.025 0.975]);
                end
                figure(29)
                grid off
                hold on
                subplot(2,3,v)
                hold on
                [ph,msg] = jbfill(Time,credible(:,1)', credible(:,2)',purp, purp,1,0.2);
                hold on
              end
            lgd = legend('model prediction', 'data', 'error bar', '95% credible interval',  'Location', 'SouthEast');
            lgd.FontSize = 14;

    end
    if GLU>0 && LPS>0

             vars = [23,24,27];

             subplot(2,3,1)
             % --- Prediction
             plot(T_best, Y_best(:,vars(1)), 'color', 'k', 'LineWidth', 1.2)
             hold on
             %plot(Time, mean_y(:,1), 'color', purp, 'LineWidth', 0.8)
             %hold on
             % --- Validation data
             scatter(timev(:,1), v_data(:,1), 30, 'ro', 'MarkerFaceColor', 'r')
             yerr = abs(ev_data(:,1) - v_data(:,1));
             errorbar(timev(:,1), v_data(:,1), yerr, 'ro');
             hold on                      
             xlabel('Time (hour)');
             ylabel(append(speciesNames([23]), ' activity'));
             ylim([0,1.05])
              xlim([0,48])
              xticks([0 12 24 36 48])
             set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Ayala et al. (2019)'); t.FontSize = 14;
             grid off
            
             subplot(2,3,2)
             % --- Prediction
             plot(T_best, Y_best(:,vars(2)), 'color', 'k', 'LineWidth', 1.2)
             hold on
             %plot(Time, mean_y(:,2), 'color', purp, 'LineWidth', 0.8)
             %hold on
             % --- Validation data
             scatter(timev(:,3), v_data(:,3),  30, 'ro', 'MarkerFaceColor', 'r')
             yerr = abs(ev_data(:,3) - v_data(:,3));
             errorbar(timev(:,3), v_data(:,3),  yerr, 'ro');
             hold on             
            
             xlabel('Time (hour)');
             ylabel(append(speciesNames([24]), ' activity'));
             ylim([0,1.05])
             xlim([0,48])
             xticks([0 12 24 36 48])
             set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Ayala et al. (2019)'); t.FontSize = 14;             
             grid off
          
         
             subplot(2,3,3);
             % --- Prediction
             plot(T_best, Y_best(:,vars(3)),'color', 'k', 'LineWidth', 1.2)
             hold on
             %plot(Time, mean_y(:,9),'color', purp, 'LineWidth', 0.8)
             %hold on
             % --- Validation data
             scatter(timev(:,4), v_data(:,4), 30, 'ro', 'MarkerFaceColor', 'r')
             yerr = abs(ev_data(:,4)-v_data(:,4));
             errorbar(timev(:,4), v_data(:,4), yerr,'ro');
             hold on
             xlabel('Time (hour)');
             ylabel(append(speciesNames([27]), ' activity')); 
             ylim([0,1.05])
             xlim([0,48])
             xticks([0 12 24 36 48])
             set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Yang et al. (2008)'); t.FontSize = 14;
             lgd = legend('model', 'model (mean)', 'CI-95%', 'data', '',  'Location', 'SouthEast');
             lgd.FontSize = 14;
             grid off

             
              for v = 1:length(vars)
                for t = 1:length(Time)
                    credible(t,:) = quantile(abs(Yp(:,t,vars(v))), [0.025 0.975]);
                end
                figure(29)
                grid off
                hold on
                subplot(2,3,v)
                hold on
                [ph,msg] = jbfill(Time,credible(:,1)', credible(:,2)',purp, purp,1,0.2);
                hold on
              end
            lgd = legend('model prediction', 'data', 'error bar', '95% credible interval',  'Location', 'SouthEast');
            lgd.FontSize = 14;
    end


end
end


 
