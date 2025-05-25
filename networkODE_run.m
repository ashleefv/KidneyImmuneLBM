function [T_best, Y_best] = networkODE_run(tspan, y0, params, tau_index, k_index, n_index, W_index, mode)

warning('off', 'all')

load invitro_data.mat

GLU = params{1}(1,1);
LPS = params{1}(1,2);

Time = [0:1:48]; % time span

speciesNames = params{4};

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

    % load posteriors of prediction
    load data/prediction_posterior_GLU.mat
end 

if GLU==0 && LPS>0                                   % condition (LPS = 1, GLU = 0)
   
    t_data = train_data([3,7,11,15,19,23],:);        % Training data for given time points
    e_data = Terror([3,7,11,15,19,23],:);            % Error bars for training
    v_data = v_data([3,7,11,15,19],:);               % Validation data for given time points
    ev_data = Verror([3,7,11,15,19],:);              % Error bars for validation
    timedata = time_data([3,7,11,15,19,23],:);       % Time points for training set
    timev = vtime([3,7,11,15,19],:);                 % Time points for validation set

    % load posteriors of prediction
    load data/prediction_posterior_LPS.mat
end

if GLU>0 && LPS>0                                    % condition (GLU = 1, LPS = 1)
    t_data = train_data([4,8,12,16,20,24],:);        % Training data for given time points
    e_data = Terror([4,8,12,16,20,24],:);            % Error bars for training
    v_data = v_data([4,8,12,16,20],:);               % Validation data for given time points
    ev_data = Verror([4,8,12,16,20],:);              % Error bars for validation
    timedata = time_data([4,8,12,16,20,24],:);       % Time points for training set
    timev = vtime([4,8,12,16,20],:);                 % Time points for validation set
    
    % load posteriors of prediction
    load data/prediction_posterior_both.mat
end


options=[];
[T_best, Y_best] = ode23s(@networkODE,tspan,y0,options,params);


%% plot mean (acceptable) parameters for each treatment condition
% fitting plots

if mode==1

            % matching figure numbers to main article
            if GLU>0 && LPS==0
                figure(4)
                filename = 'Fig4';
            elseif GLU==0 && LPS>0
                figure(5)
                filename = 'Fig5';
            else
                figure(6)
                filename = 'Fig6';
            end

            orange = [245 119 41]./256; % color
            hold on
            vars = [23, 24, 25, 6, 13, 22, 20, 12, 27, 28];

             subplot(2,4,1)
             hold on
             % --- Prediction
             plot(T_best, Y_best(:,23), 'color', 'k',  'LineWidth', 1.2)
             hold on
             % --- Training data
             scatter(timedata(:,1), t_data(:,1), 30, 'o', 'MarkerFaceColor', orange,'MarkerEdgeColor',orange)
             yerr = abs(e_data(:,1)-t_data(:,1));
             errorbar(timedata(:,1), t_data(:,1), yerr, 'o','MarkerEdgeColor',orange);
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
             scatter(timedata(:,3), t_data(:,3),  30, 'o', 'MarkerFaceColor', orange,'MarkerEdgeColor',orange)
             yerr = abs(e_data(:,3)-t_data(:,3));
             errorbar(timedata(:,3), t_data(:,3),  yerr, 'o','MarkerEdgeColor',orange);
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
             scatter(timedata(:,2), t_data(:,2),  30, 'o', 'MarkerFaceColor', orange,'MarkerEdgeColor',orange)
             yerr = abs(e_data(:,2)-t_data(:,2));
             errorbar(timedata(:,2), t_data(:,2), yerr, 'o','MarkerEdgeColor',orange);
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
             scatter(timedata(:,4), t_data(:,4), 30, 'o', 'MarkerFaceColor', orange,'MarkerEdgeColor',orange)
             yerr = abs(e_data(:,4)-t_data(:,4));
             errorbar(timedata(:,4), t_data(:,4), yerr, 'o','MarkerEdgeColor',orange);
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
             scatter(timedata(:,5), t_data(:,5),  30, 'o', 'MarkerFaceColor', orange,'MarkerEdgeColor',orange)
             yerr = abs(e_data(:,5)-t_data(:,5));
             errorbar(timedata(:,5), t_data(:,5), yerr, 'o','MarkerEdgeColor',orange);
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
             scatter(timedata(:,6), t_data(:,6),  30, 'o', 'MarkerFaceColor', orange,'MarkerEdgeColor',orange)
             yerr = abs(e_data(:,6)-t_data(:,6));
             errorbar(timedata(:,6), t_data(:,6), yerr, 'o','MarkerEdgeColor',orange);
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
              scatter(timedata(:,7),  t_data(:,7), 30, 'o', 'MarkerFaceColor', orange,'MarkerEdgeColor',orange)
              yerr = abs(e_data(:,7)-t_data(:,7));
              errorbar(timedata(:,7),  t_data(:,7), yerr, 'o','MarkerEdgeColor',orange);
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
              scatter(timedata(:,9),  t_data(:,9), 30, 'o', 'MarkerFaceColor', orange,'MarkerEdgeColor',orange)
              yerr = abs(e_data(:,9)-t_data(:,9));
              errorbar(timedata(:,9),  t_data(:,9), yerr, 'o','MarkerEdgeColor',orange);
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

              blue = 	[0 0.4470 0.7410];
              labelstring = {'a', 'b', 'c', 'd', 'e','f','g','h'};

              for v = 1:8
                  for t = 1:length(Time)
                      credible(t,:) = quantile(abs(Yp(:,t,vars(v))), [0.025 0.975]);
                  end
                  
                  hold on
                  subplot(2,4,v)
                  hold on
                  [ph,msg] = jbfill(Time,credible(:,1)', credible(:,2)',blue,blue,1,0.2);

                  text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 12)
                  hold on    

              end
              subplot(2,4,1)
              lgd = legend('model prediction', 'data', 'error bar', '95% credible interval',  'Location', 'SouthEast');
              lgd.FontSize = 12;
                
widthInches = 20;
ScriptForExportingImagesForAJP

elseif mode==2
    % validation plots

    purp = [0.4940    0.1840    0.5560]; % color
    

    labelstring = {'a', 'b', 'c', 'd', 'e','f'};

    if GLU>0 && LPS==0
             figure(7)
            filename = 'Fig7';
 labelstring = {'a', 'b', 'c', 'd', 'e','f'};
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
             plot(T_best, Y_best(:,vars(4)), 'color', 'k', 'LineWidth', 1.2)
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
             plot(T_best, Y_best(:,vars(5)), 'color', 'k', 'LineWidth', 1.2)
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
                figure(7)
                grid off
                hold on
                subplot(2,3,v)
                hold on
                [ph,msg] = jbfill(Time,credible(:,1)', credible(:,2)',purp, purp,1,0.2);
                hold on
                text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 12)
            end
            subplot(2,3,1)
            lgd = legend('model prediction', 'data', 'error bar', '95% credible interval',  'Location', 'SouthEast');
            lgd.FontSize = 14;

            
    end
    if GLU==0 && LPS>0
        figure(8)
filename = 'Fig8';
 labelstring = {'a', 'b', 'c', 'd', 'e','f'};
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
                figure(8)
                grid off
                hold on
                subplot(2,3,v)
                hold on
                [ph,msg] = jbfill(Time,credible(:,1)', credible(:,2)',purp, purp,1,0.2);
                hold on
                text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 12)
            end
            subplot(2,3,1)
            lgd = legend('model prediction', 'data', 'error bar', '95% credible interval',  'Location', 'SouthEast');
            lgd.FontSize = 14;

    end
    if GLU>0 && LPS>0
             figure(9)
             filename = 'Fig9';
              labelstring = {'a', 'b', 'c'};
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
             grid off

             
             for v = 1:length(vars)
                for t = 1:length(Time)
                    credible(t,:) = quantile(abs(Yp(:,t,vars(v))), [0.025 0.975]);
                end
                figure(9)
                grid off
                hold on
                subplot(2,3,v)
                hold on
                [ph,msg] = jbfill(Time,credible(:,1)', credible(:,2)',purp, purp,1,0.2);
                hold on
                text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 12)
            end
            subplot(2,3,1)
             lgd = legend('model prediction', 'data', 'error bar', '95% credible interval',  'Location', 'SouthEast');
             lgd.FontSize = 14;
    end
widthInches = 15;
ScriptForExportingImagesForAJP
elseif  mode == 3
% regulatory nodes

    linest = ["-", "--", "-.", ":", "--.", "-", "--", "-.", ":"];
    
    if GLU>0 && LPS==0
        figure(10)
        filename = 'Fig10';
         labelstring = {'a', 'b', 'c', 'd', 'e'};
    elseif GLU==0 && LPS>0
        figure(11)
        filename = 'Fig11';
         labelstring = {'a', 'b', 'c', 'd', 'e'};
    elseif GLU>0 && LPS>0
        figure(12)
        filename = 'Fig12';
         labelstring = {'a', 'b', 'c', 'd', 'e'};
    end
    
    var = [3,8,10];
    for j = 1:length(var)
        subplot(2,3,1)
        plot(T_best, Y_best(:,var(j)), linest(j), 'LineWidth', 2)
        hold on
        grid off
        ylim([0,1.05])
        xlim([0,48])
        xticks([0 12 24 36 48])
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
        grid off
        ylim([0,1.05])
        xlim([0,48])
        xticks([0 12 24 36 48])
        ylabel("Activity"); xlabel('Time (hour)')
        set(gca,'FontName','Arial','FontSize',18);
        lgd = legend(params{4}([9,14,15,19]), 'Location', 'SouthEast');
        lgd.FontSize = 14;
        
    end
    var = [7,11,12, 18];
    for j=1:length(var)
        subplot(2,3,3)
        plot(T_best, Y_best(:,var(j)), linest(j), 'LineWidth', 2); 
        hold on
        grid off
        ylim([0,1.05])
        xlim([0,48])
        xticks([0 12 24 36 48])
        ylabel("Activity"); xlabel('Time (hour)')
        set(gca,'FontName','Arial','FontSize',18);
        lgd = legend(params{4}([7,11,12,18]), 'Location', 'SouthEast');
        lgd.FontSize = 14;
    end
    var = [4,5,26,16,17];
    for j=1:length(var)
        subplot(2,3,4)
        plot(T_best, Y_best(:,var(j)), linest(j), 'LineWidth', 2); 
        hold on
        grid off
        ylim([0,1.05])
        xlim([0,48])
        xticks([0 12 24 36 48])
       
        ylabel("Activity"); xlabel('Time (hour)')
        set(gca,'FontName','Arial','FontSize',18);
        lgd = legend(params{4}([4,5,26,16,17]), 'Location', 'SouthEast');
        lgd.FontSize = 14;
    end
    var = [28,29];
    for j=1:length(var)
        subplot(2,3,5)
        plot(T_best, Y_best(:,var(j)), linest(j), 'LineWidth', 2)
        hold on
        grid off
        ylim([0,1.05])
        xlim([0,48])
        xticks([0 12 24 36 48])
        ylabel("Activity"); xlabel('Time (hour)')
        set(gca,'FontName','Arial','FontSize',18);
        lgd = legend(params{4}([28,29]), 'Location', 'SouthEast');
        lgd.FontSize = 14;
    end

for v = 1:5
    subplot(2,3,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 12)
end
widthInches = 15;
ScriptForExportingImagesForAJP
end