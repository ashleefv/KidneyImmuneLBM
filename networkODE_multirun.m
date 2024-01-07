function networkODE_multirun(tspan, y0, params, p_fitted, global_p_best, tau_index, k_index, n_index, W_index)

warning('off', 'all')

load invitro_data.mat

GLU = params{1}(1,1);
LPS = params{1}(1,2);
Time = [0:1:48]; % time span

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

purp = [0.4940    0.1840    0.5560];

speciesNames = params{4};


options=[];
[T_best, Y_best] = ode23s(@networkODE,tspan,y0,options,params);


%%
dp = params;

dp{2}(tau_index) = global_p_best(1:size_tau);
dp{1}(1,W_index) = global_p_best(size_tau+1:size_tau+size_W);
dp{1}(2,n_index) = global_p_best(size_tau+size_W+1:size_tau+size_W+size_n);
dp{1}(3,k_index) = global_p_best(size_tau+size_W+size_n+1:size_tau+size_W+size_n+size_k);

%% plot global best-fit parameters for each treatment condition


if mode==1

figure(13)
hold on

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
             ylabel(speciesNames([23]));
             ylim([0,1.5]);
             set(gca,'FontName','Arial','FontSize',18);
             t = title('Tessaro et al. (2020)'); t.FontSize = 14;
             grid on
                       
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
             xlabel('Time (hour)');
             ylabel(speciesNames([24]));
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
             t = title('Tessaro et al. (2020)'); t.FontSize = 14;
             grid on
          
     
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
             xlabel('Time (hour)');
             ylabel(speciesNames([25])); 
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18); 
             t = title('Li et al. (2021)'); t.FontSize = 14;
             grid on

         
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
             xlabel('Time (hour)');
             ylabel(speciesNames([6])); 
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18); 
             t = title('Li et al. (2021)'); t.FontSize = 14;
             grid on

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
             xlabel('Time (hour)');
             ylabel(speciesNames([13])); 
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
             t = title('de Souza et al. (2007)'); t.FontSize = 14;
             grid on

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
             xlabel('Time (hour)');
             ylabel(speciesNames([22])); 
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
             t = title('Nakagawa et al. (2006)'); t.FontSize = 14;

             grid on

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
              xlabel('Time (hour)');
              ylabel(speciesNames([20])); 
              ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
              t = title('Kasai et al. (2005)'); t.FontSize = 14;
              grid on

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
              xlabel('Time (hour)');
              ylabel(speciesNames([12]));
              ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);  
              hold on
              t = title('Yang et al. (2008)'); t.FontSize = 14;

             
            
              grid on
        
end

%%
if mode==2
    vars = [23, 24, 25, 6, 13, 22, 20, 12, 27, 28];

    figure(14)

    if GLU>0 && LPS==0

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
             ylabel(speciesNames([23]));
             ylim([0,1.5]);
             set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Ayala et al. (2019)'); t.FontSize = 14;
             grid on
            
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
             ylabel(speciesNames([24]));
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Ayala et al. (2019)'); t.FontSize = 14;             
             grid on
          
         
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
             ylabel(speciesNames([25])); 
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('He et al. (2020)'); t.FontSize = 14;
             grid on

         
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
             ylabel(speciesNames([27])); 
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Yang et al. (2008)'); t.FontSize = 14;
             grid on


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
              ylabel(speciesNames([12]));
              ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);  
              hold on
              t = title('Seo et al. (2022)'); t.FontSize = 14;
              hold on
            
              grid on

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
              ylabel(speciesNames([20])); 
              ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
              hold on
              t = title('Ramanathan et al. (2002)'); t.FontSize = 14;
              
              grid on

            vars = [23,24,25,27,12,20];
            for v = 1:length(vars)
                for t = 1:length(Time)
                    credible(t,:) = quantile(abs(Yp(:,t,vars(v))), [0.025 0.975]);
                end
                figure(14)
                grid on
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
             ylabel(speciesNames([23]));
             ylim([0,1.5]);
             set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Ayala et al. (2019)'); t.FontSize = 14;
             grid on
            
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
             ylabel(speciesNames([24]));
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Ayala et al. (2019)'); t.FontSize = 14;             
             grid on
          
         
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
             ylabel(speciesNames([25])); 
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('He et al. (2020)'); t.FontSize = 14;
             grid on

         
             subplot(2,3,4);
             % --- Prediction
             plot(T_best, Y_best(:,vars(9)),'color', 'k', 'LineWidth', 1.2)
             hold on
             %plot(Time, mean_y(:,9),'color', purp, 'LineWidth', 0.8)
             %hold on
             % --- Validation data
             scatter(timev(:,4), v_data(:,4), 30, 'ro', 'MarkerFaceColor', 'r')
             yerr = abs(ev_data(:,4)-v_data(:,4));
             errorbar(timev(:,4), v_data(:,4), yerr,'ro');
             hold on
             xlabel('Time (hour)');
             ylabel(speciesNames([27])); 
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Yang et al. (2008)'); t.FontSize = 14;
             grid on


             subplot(2,3,5)
             % --- Prediction
              plot(T_best, Y_best(:,vars(10)), 'color', 'k', 'LineWidth', 1.2)
              hold on
              %plot(Time, mean_y(:,10), 'color', purp, 'LineWidth', 0.8)
              %hold on   
              % --- Validation data
              scatter(timev(:,8), v_data(:,8), 30, 'ro', 'MarkerFaceColor', 'r')
              yerr = abs(ev_data(:,8)-v_data(:,8));
              errorbar(timev(:,8), v_data(:,8), yerr, 'ro');
              hold on              
              xlabel('Time (hour)');
              ylabel(speciesNames([28]));
              ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);  
              hold on
              t = title('Seo et al. (2022)'); t.FontSize = 14;
              hold on
            
              grid on

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
              ylabel(speciesNames([20])); 
              ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
              hold on
              t = title('Ramanathan et al. (2002)'); t.FontSize = 14;
              grid on

              vars = [23,24,25,27,28,20];
              for v = 1:length(vars)
                for t = 1:length(Time)
                    credible(t,:) = quantile(abs(Yp(:,t,vars(v))), [0.025 0.975]);
                end
                figure(14)
                grid on
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
             ylabel(speciesNames([23]));
             ylim([0,1.5]);
             set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Ayala et al. (2019)'); t.FontSize = 14;
             grid on
            
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
             ylabel(speciesNames([24]));
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Ayala et al. (2019)'); t.FontSize = 14;             
             grid on
          
         
             subplot(2,3,3);
             % --- Prediction
             plot(T_best, Y_best(:,vars(9)),'color', 'k', 'LineWidth', 1.2)
             hold on
             %plot(Time, mean_y(:,9),'color', purp, 'LineWidth', 0.8)
             %hold on
             % --- Validation data
             scatter(timev(:,4), v_data(:,4), 30, 'ro', 'MarkerFaceColor', 'r')
             yerr = abs(ev_data(:,4)-v_data(:,4));
             errorbar(timev(:,4), v_data(:,4), yerr,'ro');
             hold on
             xlabel('Time (hour)');
             ylabel(speciesNames([27])); 
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title(''); t.FontSize = 14;
             lgd = legend('model', 'model (mean)', 'CI-95%', 'data', '',  'Location', 'SouthEast');
             lgd.FontSize = 14;
             grid on

             vars = [23,24,27];
              for v = 1:length(vars)
                for t = 1:length(Time)
                    credible(t,:) = quantile(abs(Yp(:,t,vars(v))), [0.025 0.975]);
                end
                figure(14)
                grid on
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
%% regulatory nodes
% Plot regulatory nodes in the network to study the dynamics in absence of
% training or validation data


if mode == 3
    linest = ["-", "--", "-.", ":", "--."];
    figure(10);
    
    var = [3,8,10];
    for j = 1:length(var)
        subplot(2,3,1)
        plot(T_best, Y_best(:,var(j)), linest(j), 'LineWidth', 2)
        hold on
        grid on
        ylabel("Activity"); xlabel('Time (hour)')
        set(gca,'FontName','Arial','FontSize',18);
        hold on
     lgd = legend(params{4}([3,8,10]), 'Location', 'SouthEast');
     lgd.FontSize = 14;

    end
    var = [9,14,15,19];
    for j = 1:length(var)
        subplot(2,3,2)
        plot(T_best, Y_best(:,var(j)), linest(j), 'LineWidth', 2); 
        hold on
        grid on
        ylabel("Activity"); xlabel('Time (hour)')
        set(gca,'FontName','Arial','FontSize',18);
        lgd = legend(params{4}([9,14,15,19]), 'Location', 'SouthEast');
        lgd.FontSize = 14;
        
    end
    var = [7,11,18];
    for j=1:length(var)
        subplot(2,3,3)
        plot(T_best, Y_best(:,var(j)), linest(j), 'LineWidth', 2); 
        hold on
        grid on
        ylabel("Activity"); xlabel('Time (hour)')
        set(gca,'FontName','Arial','FontSize',18);
        lgd = legend(params{4}([7,11,18]), 'Location', 'SouthEast');
        lgd.FontSize = 14;
    end
    var = [4,5,26,16,17];
    for j=1:length(var)
        subplot(2,3,4)
        plot(T_best, Y_best(:,var(j)), linest(j), 'LineWidth', 2); 
        hold on
        grid on
        ylabel("Activity"); xlabel('Time (hour)')
        set(gca,'FontName','Arial','FontSize',18);
        lgd = legend(params{4}([4,5,26,16,17]), 'Location', 'SouthEast');
        lgd.FontSize = 14;
    end
    var = [28,29,30];
    for j=1:length(var)
        subplot(2,3,5)
        plot(T_best, Y_best(:,var(j)), linest(j), 'LineWidth', 2)
        hold on
        grid on
        ylabel("Activity"); xlabel('Time (hour)')
        set(gca,'FontName','Arial','FontSize',18);
        lgd = legend(params{4}([28,29,30]), 'Location', 'SouthEast');
        lgd.FontSize = 14;
    end

end
%% Predictions on trained model


linest = ["-", "--", "-.", ":", "--."];
R1  = [0,0.25,0.5,0.75,1];
for i = [1:length(R1)]
    
     
    %params{1}(1,[9,11,25,34]) = R1(i);
    
   
    params{1}(1,1) = R1(i);
    

    %params{1}(1,37) = R1(i);
    %params{3}(29) = YM2(i);
    %params{1}(12) = YM1(i);

    
    %figure(i); plot(Tout, Yout(:,[20, 24, 25, 26, 27, 30])); 
   
    %legendInfo{i} = ['W(ROS) = ' num2str(R1(i))];
    
    %VEGF_mean(i) = real(mean(Yout(:,27)));
    %IL6_mean(i) = real(mean(Yout(:,23)));
    %NO_mean(i) = real(mean(Yout(:,20)));
    %gap_end(i) = real((Yout(end,30)));
    %ROS_mean(i) = real(mean(Yout(:,12)));
    %Ca_mean(i) = real(mean(Yout(:,29)));
    %junc_mean(i) = real(mean(Yout(:,28)));
    
    options = [];
    [To, Yo] = ode23s(@networkODE,tspan,y0,options,params);

    figure(21)
    subplot(2,2,1)
    plot(To, Yo(:,23),  linest(i), 'LineWidth', 2)
    hold on
    ylabel("IL-6 activity"); xlabel('Time(hour)')
    set(gca,'FontName','Arial','FontSize',18);
    grid on;
    hold on

    subplot(2,2,2)
    plot(To, Yo(:,13),  linest(i), 'LineWidth', 2)
    ylabel("ROS_e_c activity"); xlabel('Time(hour)')
    set(gca,'FontName','Arial','FontSize',18);
    grid on;
    hold on
    
    subplot(2,2,3)
    plot(To, Yo(:,26),  linest(i), 'LineWidth', 2)
    hold on
    ylabel("PLC-\gamma activity"); xlabel('Time(hour)')
    set(gca,'FontName','Arial','FontSize',18);
    grid on;
    hold on

    subplot(2,2,4)
    plot(To, Yo(:,30),  linest(i), 'LineWidth', 2)
    ylabel("Change in Gap Width"); xlabel('Time(hour)')
    set(gca,'FontName','Arial','FontSize',18);
    grid on;
    hold on
    legendI{i} = ['W_G_L_U = ' num2str(R1(i))];
    legend(legendI)
    hold on

    params{1}(1,1) = 1;
    params{1}(1,2) = R1(i);
    
    figure(22)
    subplot(2,2,1)
    plot(To, Yo(:,24),  linest(i), 'LineWidth', 2)
    hold on
    ylabel("TNF-\alpha activity"); xlabel('Time(hour)')
    set(gca,'FontName','Arial','FontSize',18);
    grid on;
    hold on

    subplot(2,2,2)
    plot(To, Yo(:,25),  linest(i), 'LineWidth', 2)
    ylabel("IL-1\beta activity"); xlabel('Time(hour)')
    set(gca,'FontName','Arial','FontSize',18);
    grid on;
    hold on
    
    subplot(2,2,3)
    plot(To, Yo(:,26),  linest(i), 'LineWidth', 2)
    hold on
    ylabel("PLC-\gamma activity"); xlabel('Time(hour)')
    set(gca,'FontName','Arial','FontSize',18);
    grid on;
    hold on

    subplot(2,2,4)
    plot(To, Yo(:,20),  linest(i), 'LineWidth', 2)
    ylabel("NO activity"); xlabel('Time(hour)')
    set(gca,'FontName','Arial','FontSize',18);
    grid on;
    hold on
    legendInfo{i} = ['W_G_L_U, W_L_P_S=(1,' num2str(R1(i)) ')'];
    hold on
    

   

end

hold on
legend(legendInfo)


%% plot ensemble predictions for all fitted parameter sets (100)
if mode==4

% Uncomment next two lines and read fitted params for each treatment condition (#: GLU, LPS, or both) 
% p_fitted = readmatrix('~/#_fitted_params.csv'); p_fitted = p_fitted([2:end],[1:end-1]);
% error_fitted = readmatrix('~/#_fitted_params.csv'); error_fitted = error_fitted([2:end],end);

col = [0 0.4470 0.7410];
Params_error = [error_fitted p_fitted]; % coef_fitted is the parameter output from all the LHS fmincon runs
error_column = 1; % first column stores error_fitted
[sorted_Params_error, index] = sortrows(Params_error,error_column);

param_posterior = sorted_Params_error(:,[2:end]);

figure(13)

for i = 1:100
            dp = params;
            dp{2}(tau_index) = param_posterior(i,[1:size_tau]);
            %dp{1}(1,W_index) = param_posterior(i,[size_tau+1:size_tau+size_W]);
            %dp{1}(2,n_index) = param_posterior(i,[size_tau+size_W+1:size_tau+size_W+size_n]);
            %dp{1}(3,k_index) = param_posterior(i,[size_tau+size_W+size_n+1:end]);

            options = [];
            % options = odeset('RelTol',1e-20); %,'MaxStep',1e-16);
            [Tout, Yout] = ode23s(@networkODE,tspan,y0,options,dp);
    
             subplot(2,4,1)
             hold on
             plot(Tout, Yout(:,23), '--', 'color', col)
             hold on
                       
             subplot(2,4,2)
             hold on
             % --- Prediction
             plot(Tout, Yout(:,24), '--', 'color', col)
             hold on
             
          
     
             subplot(2,4,3)
             hold on
             % --- Prediction
             plot(Tout, Yout(:,25), '--', 'color', col)
             hold on
             

         
             subplot(2,4,4)
             hold on
             % --- Prediction
             plot(Tout, Yout(:,6), '--', 'color', col)
             hold on
             
             subplot(2,4,5);
             % --- Prediction
             hold on
             plot(Tout, Yout(:,13), '--', 'color', col)
             hold on
             

             subplot(2,4,6);
             % --- Prediction
             hold on
             plot(Tout, Yout(:,22), '--', 'color', col)
             hold on
             

             subplot(2,4,7)
             hold on
             % --- Prediction
             plot(Tout, Yout(:,20), '--', 'color', col)
             hold on
             

             subplot(2,4,8);
             hold on
             % --- Prediction
             plot(Tout, Yout(:,12), '--', 'color', col)
             hold on
             
end
%% Monte Carlo Simulation
Params_error = [error_fitted p_fitted]; % coef_fitted is the parameter output from all the LHS fmincon runs
error_column = 1; % first column stores error_fitted
[sorted_Params_error, index] = sortrows(Params_error,error_column);



figure(20)
%-- generate figure to show the error of all the parameter values from LHS
for i = 1:size(p_fitted, 2)
    plot(1:100,sorted_Params_error(:,error_column),'o')
    hold on
end


%%
dp = params;

p_posterior = sorted_Params_error([1:100],[2:end]);
Ns = 5000;
population = zeros(1,length(global_p_best));
%population = zeros(1,size(p_fitted,2));
s = RandStream('mlfg6331_64');
tic
for j = 1:Ns
    for m = 1:length(global_p_best)
        population(j,m) = randsample(s, p_posterior(:,m), 1);
    end
end
%%
for j = 1:Ns
    dp{2}(tau_index) = population(j,(1:size_tau));
    %dp{1}(1,W_index) = population(j,(size_tau+1:size_tau+size_W));
    %dp{1}(2,n_index) = population(j,(size_tau+size_W+1:size_tau+size_W+size_n));
    %dp{1}(3,k_index) = population(j,(size_tau+size_W+size_n+1:end)); 
    
    options = [];
    % options = odeset('RelKTol',1e-20); %,'MaxStep',1e-16);
    [Tout, Yout] = ode23s(@networkODE,tspan,y0,options,dp);
    Yp(j,:,:) = Yout;
    dp = params;
   
end
toc
%%
Time = [0:1:48];
vars = [23, 24, 25, 6, 13, 22, 20, 12, 27, 28];

for v = 1:length(vars)
    for i = 1:length(Time)
    
        mean_y(i,v) = mean(abs(Yp(:,i,vars(v))));
        sd_y(i,v) = std(abs(Yp(:,i,vars(v))));
        

    end
    ts = tinv([0.025  0.975],length(Time)-1);
    CI95 = abs(mean_y(:,v) + ts.*(sd_y(:,v)/sqrt(length(Time))));
    lower_CI(:,v) = CI95(:,1);
    upper_CI(:,v) = CI95(:,2);

    
end

%% plot mean posterior predictions and 95% confidence intervals
         
% save('workspaces/both_best_fitted_MC')
%%

for v = 1:8
    for t = 1:length(Time)
    
       credible(t,:) = quantile(abs(Yp(:,t,vars(v))), [0.025 0.975]);
    end
    figure(13)
    grid on
    hold on
    subplot(2,4,v)
    hold on
    % --- Prediction
    %plot(Time, mean_y(:,v), 'color', 'b', 'LineWidth', 0.8)
    %hold on
    [ph,msg] = jbfill(Time,credible(:,1)', credible(:,2)','b','b',1,0.2);
    hold on
    ylabel(params{4}(vars(v)))
    xlabel('Time (hr)')
    

end

lgd = legend('model prediction', 'data', 'error bar', '95% credible interval',  'Location', 'SouthEast');
lgd.FontSize = 12;

end

end

 
