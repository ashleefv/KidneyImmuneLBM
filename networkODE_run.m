function [Tout, Yout] = networkODE_run(tspan, y0, params, mode)

load invitro_data.mat

options = [];
% opt = odeset('RelTol',1e-20); %,'MaxStep',1e-16);

[t, y] = ode23s(@networkODE,tspan,y0,options,params);


Yout = real(y);
Tout = t;


GLU = params{1}(1,1);
LPS = params{1}(1,2);
Time = [0:1:48]; % time span



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
%%
purp = [0.4940    0.1840    0.5560];

if mode == 1
    Y_var = Yout(:,[23,24,25,6,13,22,20,12]);        % Training variables
    num_outputvar = size(Y_var,2);                   % number of variables
    speciesNames = params{4};                         
    figure(4)
     for i = 1:num_outputvar
         for t=1:length(Time)
             p = kde(Y_var(t,i), 0.05);
             meanholder(t) = meanKDE(p);             % returns mean of KDE
             sdholder(t) = sqrt(covar(p));           % returns variance of KDE
             q = quantile(sample(p,10000),[0.05 0.25 0.75 0.95]); % when sample() is a vector, q(i) is the p(i)th quartiles (0.05 or 0.25 or ...), q is same size as columns of sample()
             qholder(t,i,:) = q;
         end
         % convert any possible negative quantile values to a very small number, 1e-100
         for k = 1:4                                 % determined by the confidences you want to calculate (1:2 if you only want one CI say 90th with quantile 0.05 and 0.95)
             for iS = 1:length(Time)
                 if qholder(iS,i,k) < 0
                     qholderPositive(iS,i,k) = 1e-100;
                 else
                     qholderPositive(iS,i,k) = qholder(iS,i,k);
                 end
             end
         end
         
         if i == 1
             col = [0 0.4470 0.7410]; 
             subplot(2,4,i)
             % --- Prediction
             plot(Tout, Y_var(:,i), 'LineWidth', 2, 'color', col)
             hold on
             % --- Training data
             scatter(timedata(:,1), t_data(:,1), 50, 'filled')
             yerr = abs(e_data(:,1)-t_data(:,1));
             errorbar(timedata(:,1), t_data(:,1), yerr, 'o');
             hold on                       
             % --- shade 50% CI
             % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',col,col,0,0.2); 
             hold on
             %--- shade 90% CI
             [ph,msg]=jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',col,col,1,0.2);
             hold on            
             xlabel('Time (hour)');
             ylabel(speciesNames([23]));
             ylim([0,1.5]);
             set(gca,'FontName','Arial','FontSize',18);
             t = title('Tessaro et al. (2020)'); t.FontSize = 14;
             grid on
             
         elseif i == 2
             col = [0 0.4470 0.7410];              
             subplot(2,4,i)
             % --- Prediction
             plot(Tout, Y_var(:,i), 'LineWidth',2, 'color', col)
             hold on
             % --- Training data
             scatter(timedata(:,3), t_data(:,3),  50, 'filled')
             yerr = abs(e_data(:,3)-t_data(:,3));
             errorbar(timedata(:,3), t_data(:,3),  yerr, 'o');
             hold on
             % --- shade 50% CI
             % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',col,col,0,0.2); 
             hold on
             % --- shade 90% CI
             [ph,msg]=jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',col,col,1,0.2);
             hold on
             xlabel('Time (hour)');
             ylabel(speciesNames([24]));
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
             t = title('Tessaro et al. (2020)'); t.FontSize = 14;
             grid on
          
         elseif i==3
             col = [0 0.4470 0.7410]; 
             subplot(2,4,i);
             % --- Prediction
             plot(Tout, Y_var(:,i), 'LineWidth', 2, 'color', col)
             hold on
             % --- Training data
             scatter(timedata(:,2), t_data(:,2),  50, 'filled')
             yerr = abs(e_data(:,2)-t_data(:,2));
             errorbar(timedata(:,2), t_data(:,2), yerr, 'o');
             hold on
             % --- shade 50% CI
             % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',col,col,0,0.2); 
             hold on
             % --- shade 90% CI
             [ph,msg]=jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',col,col,1,0.2);
             hold on
             xlabel('Time (hour)');
             ylabel(speciesNames([25])); 
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18); 
             t = title('Li et al. (2021)'); t.FontSize = 14;
             grid on

         elseif i==4
             col = [0 0.4470 0.7410];
             subplot(2,4,i);
             % --- Prediction
             plot(Tout, Y_var(:,i), 'LineWidth', 2, 'color', col)
             hold on
             % --- Training data
             scatter(timedata(:,4), t_data(:,4), 50, 'filled')
             yerr = abs(e_data(:,4)-t_data(:,4));
             errorbar(timedata(:,4), t_data(:,4), yerr, 'o');
             hold on
             % --- shade 50% CI
             % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',col,col,0,0.2);
             % --- shade 90% CI
             [ph,msg]=jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',col,col,1,0.2);
             hold on
             xlabel('Time (hour)');
             ylabel(speciesNames([6])); 
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18); 
             t = title('Li et al. (2021)'); t.FontSize = 14;
             grid on

         elseif i==5
             col = [0 0.4470 0.7410];
             subplot(2,4,i);
             % --- Prediction
             plot(Tout, Y_var(:,i), 'LineWidth', 2, 'color', col)
             hold on
             % --- Training data
             scatter(timedata(:,5), t_data(:,5),  50, 'filled')
             yerr = abs(e_data(:,5)-t_data(:,5));
             errorbar(timedata(:,5), t_data(:,5), yerr, 'o');
             hold on           
             % --- shade 50% CI
             % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',col,col,0,0.2); 
             % --- shade 90% CI
             [ph,msg]=jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',col,col,1,0.2);
             hold on
             xlabel('Time (hour)');
             ylabel(speciesNames([13])); 
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
             t = title('de Souza et al. (2007)'); t.FontSize = 14;
             grid on

         elseif i==6
             col = [0 0.4470 0.7410]; 
             subplot(2,4,i)
             % --- Prediction
             plot(Tout,Y_var(:,i), 'LineWidth', 2, 'color', col);
             hold on
             % --- Training data
             scatter(timedata(:,6), t_data(:,6),  50, 'filled')
             yerr = abs(e_data(:,6)-t_data(:,6));
             errorbar(timedata(:,6), t_data(:,6), yerr, 'o');
             hold on             
             % --- shade 50% CI
             % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',col,col,0,0.2);
             % --- shade 90% CI
             [ph,msg]=jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',col,col,1,0.2);
             hold on
             xlabel('Time (hour)');
             ylabel(speciesNames([22])); 
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
             t = title('Nakagawa et al. (2006)'); t.FontSize = 14;

             grid on

         elseif i==7
              col = [0 0.4470 0.7410];
              subplot(2,4,i)
              % --- Prediction
              plot(Tout,Y_var(:,i), 'LineWidth', 2, 'color', col);
              hold on
              % --- Training data
              scatter(timedata(:,7),  t_data(:,7), 50, 'filled')
              yerr = abs(e_data(:,7)-t_data(:,7));
              errorbar(timedata(:,7),  t_data(:,7), yerr, 'o');
              hold on     
              % --- shade 50% CI
              % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',col,col,0,0.2); 
              hold on
              % --- shade 90% CI
              [ph,msg]=jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',col,col,1,0.2);
              hold on
              xlabel('Time (hour)');
              ylabel(speciesNames([20])); 
              ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
              t = title('Kasai et al. (2005)'); t.FontSize = 14;
              grid on

         else
              col = [0 0.4470 0.7410];
              subplot(2,4,i)
              % --- Prediction
              plot(Tout,Y_var(:,i), 'LineWidth', 2, 'color', col);
              hold on
              % --- Training data
              scatter(timedata(:,9),  t_data(:,9), 50, 'filled')
              yerr = abs(e_data(:,9)-t_data(:,9));
              errorbar(timedata(:,9),  t_data(:,9), yerr, 'o');
              hold on
              % --- shade 50% CI
              % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',col,col,0,0.2); 
              % --- shade 90% CI
              [ph,msg] = jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',col,col,1,0.2);
              hold on              
              xlabel('Time (hour)');
              ylabel(speciesNames([12]));
              ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);  
              hold on
              t = title('Yang et al. (2008)'); t.FontSize = 14;
              lgd = legend('model', 'data', '', 'CI-90%', 'Location', 'SouthEast');
              lgd.FontSize = 14;
            
              grid on
           
  
         end
 
     end

end
%% Validation plots 
% Validation plots for each treatment condition: GLU only, LPS only, GLU and LPS

if mode==2
    if GLU==1 && LPS==0

    Y_var = Yout(:,[23,24,25,27,12,20]);             % Training variables
    num_outputvar = size(Y_var,2);                   % number of variables
    speciesNames = params{4};                         
    figure(5)
     for i = 1:num_outputvar
         for t=1:length(Time)
             p = kde(Y_var(t,i), 0.05);
             meanholder(t) = meanKDE(p);             % returns mean of KDE
             sdholder(t) = sqrt(covar(p));           % returns variance of KDE
             q = quantile(sample(p,10000),[0.05 0.25 0.75 0.95]); % when sample() is a vector, q(i) is the p(i)th quartiles (0.05 or 0.25 or ...), q is same size as columns of sample()
             qholder(t,i,:) = q;
         end
         % convert any possible negative quantile values to a very small number, 1e-100
         for k = 1:4                                 % determined by the confidences you want to calculate (1:2 if you only want one CI say 90th with quantile 0.05 and 0.95)
             for iS = 1:length(Time)
                 if qholder(iS,i,k) < 0
                     qholderPositive(iS,i,k) = 1e-100;
                 else
                     qholderPositive(iS,i,k) = qholder(iS,i,k);
                 end
             end
         end
         
         if i == 1
             subplot(2,3,i)
             % --- Prediction
             plot(Tout, Y_var(:,i), 'LineWidth', 2, 'color', purp)
             hold on
             % --- Validation data
             scatter(timev(:,1), v_data(:,1), 100, '*')
             yerr = abs(ev_data(:,1) - v_data(:,1));
             errorbar(timev(:,1), v_data(:,1), yerr, '*');
             hold on            
             % --- shade 50% CI
             % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',purp,purp,0,0.2); 
             hold on
             %--- shade 90% CI
             [ph,msg]=jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',purp,purp,1,0.2);
             hold on            
             xlabel('Time (hour)');
             ylabel(speciesNames([23]));
             ylim([0,1.5]);
             set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Ayala et al. (2019)'); t.FontSize = 14;
             grid on
             
         elseif i == 2             
             subplot(2,3,i)
             % --- Prediction
             plot(Tout, Y_var(:,i), 'LineWidth',2, 'color', purp)
             hold on
             % --- Validation data
             scatter(timev(:,3), v_data(:,3),  100, '*')
             yerr = abs(ev_data(:,3) - v_data(:,3));
             errorbar(timev(:,3), v_data(:,3),  yerr, '*');
             hold on             
             % --- shade 50% CI
             % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',purp,purp,0,0.2); 
             hold on
             % --- shade 90% CI
             [ph,msg]=jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',purp,purp,1,0.2);
             hold on
             xlabel('Time (hour)');
             ylabel(speciesNames([24]));
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Ayala et al. (2019)'); t.FontSize = 14;             
             grid on
          
         elseif i==3 
             subplot(2,3,i);
             % --- Prediction
             plot(Tout, Y_var(:,i), 'LineWidth', 2, 'color', purp)
             hold on
             % --- Validation data
             scatter(timev(:,2), v_data(:,2),  100, '*')
             yerr = abs(ev_data(:,2)-v_data(:,2));
             errorbar(timev(:,2), v_data(:,2), yerr, '*');
             hold on             
             % --- shade 50% CI
             % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',purp,purp,0,0.2); 
             hold on
             % --- shade 90% CI
             [ph,msg]=jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',purp,purp,1,0.2);
             hold on
             xlabel('Time (hour)');
             ylabel(speciesNames([25])); 
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('He et al. (2020)'); t.FontSize = 14;
             grid on

         elseif i==4
             subplot(2,3,i);
             % --- Prediction
             plot(Tout, Y_var(:,i), 'LineWidth', 2, 'color', purp)
             hold on
             % --- Validation data
             scatter(timev(:,4), v_data(:,4), 100, '*')
             yerr = abs(ev_data(:,4)-v_data(:,4));
             errorbar(timev(:,4), v_data(:,4), yerr,'*');
             hold on
             % --- shade 50% CI
             % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',purp,purp,0,0.2);
             % --- shade 90% CI
             [ph,msg]=jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',purp,purp,1,0.2);
             hold on
             xlabel('Time (hour)');
             ylabel(speciesNames([27])); 
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Yang et al. (2008)'); t.FontSize = 14;
             grid on

         elseif i==5
              subplot(2,3,i)
              % --- Prediction
              plot(Tout,Y_var(:,i), 'LineWidth', 2, 'color', purp);
              hold on
              % --- Validation data
              scatter(timev(:,5), v_data(:,5), 100, '*')
              yerr = abs(ev_data(:,5)-v_data(:,5));
              errorbar(timev(:,5), v_data(:,5), yerr, '*');
              % --- shade 50% CI
              % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',purp,purp,0,0.2); 
              % --- shade 90% CI
              [ph,msg] = jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',purp,purp,1,0.2);
              hold on              
              xlabel('Time (hour)');
              ylabel(speciesNames([12]));
              ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);  
              hold on
              t = title('Seo et al. (2022)'); t.FontSize = 14;
              hold on
            
              grid on

         elseif i==6             
              subplot(2,3,i)
              % --- Prediction
              plot(Tout,Y_var(:,i), 'LineWidth', 2, 'color', purp);
              hold on
              % --- Validation data
              scatter(timev(:,7),  v_data(:,7), 100, '*')
              yerr = abs(ev_data(:,7)-v_data(:,7));
              errorbar(timev(:,7),  v_data(:,7), yerr, '*');
              hold on              
              % --- shade 50% CI
              % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',purp,purp,0,0.2); 
              hold on
              % --- shade 90% CI
              [ph,msg]=jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',purp,purp,1,0.2);
              hold on
              xlabel('Time (hour)');
              ylabel(speciesNames([20])); 
              ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
              hold on
              t = title('Ramanathan et al. (2002)'); t.FontSize = 14;
              lgd = legend('model', 'data', '', 'CI-90%', 'Location', 'SouthEast');
              lgd.FontSize = 14;
              grid on
         end
          
     end
 % Additional Validation plots: p-Junc, Calcium, Gap for GLU=1 treatment.

 Y_var = Yout(:,[28, 29, 30]);                 % Selected variables

    num_outputvar = size(Y_var,2);             % no. of variables
    speciesNames = params{4};
     for i = 1:num_outputvar
         for t=1:length(Time)
             p = kde(Y_var(t,i), 0.05);
             meanholder(t) = meanKDE(p); 
             sdholder(t) = sqrt(covar(p)); 
             q = quantile(sample(p,10000),[0.05 0.25 0.75 0.95]); 
             qholder(t,i,:)=q;
         end
         % convert any possible negative quantile values to a very small number, 1e-100
         for k = 1:4 
             for iS = 1:length(Time)
                 if qholder(iS,i,k) < 0
                     qholderPositive(iS,i,k) = 1e-100;
                 else
                     qholderPositive(iS,i,k) = qholder(iS,i,k);
                 end
             end
         end

         figure(6)
         
          if i == 1
             col = [0.4940 0.1840 0.5560]; 
             subplot(2,3,i)
             % --- Prediction
             plot(Tout, Y_var(:,i), 'LineWidth', 2, 'color', purp)
             hold on
             % --- Validation data
             scatter(timev(:,8), v_data(:,8), 100, '*')
             yerr = abs(ev_data(:,8)-v_data(:,8));
             errorbar(timev(:,8), v_data(:,8), yerr, '*');
             hold on             
             % --- shade 50% CI
             % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',purp,purp,0,0.2); 
             % --- shade 90% CI
             [ph,msg] = jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',purp,purp,1,0.2);
             hold on
             xlabel('Time (hour)');
             ylabel(speciesNames([28]));
             ylim([0,1.5]);
             set(gca,'FontName','Arial','FontSize',18);
             grid on

         elseif i == 2
             col = [0.4940 0.1840 0.5560]; 
             subplot(2,3,i)
             % --- Prediction
             plot(Tout, Y_var(:,i), 'LineWidth',2, 'color', purp)
             hold on
             % --- Validation data
             scatter(timev(:,6), v_data(:,6),  100, '*')
             yerr = abs(ev_data(:,6)-v_data(:,6));
             errorbar(timev(:,6), v_data(:,6),  yerr, '*');
             hold on             
             % --- shade 50% CI
             % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',purp,purp,0,0.2); 
             % --- shade 90% CI
             [ph,msg] = jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',purp,purp,1,0.2);
             hold on
             xlabel('Time (hour)');
             ylabel(speciesNames([29]));
             ylim([0,1.5]);
             set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Seo et al. (2022)'); t.FontSize = 14;
             grid on
          
          else
             subplot(2,3,i);
             % --- Prediction
             plot(Tout, Y_var(:,i), 'LineWidth', 2, 'color', purp)
             hold on
             % --- Validation data
             scatter(timedata(:,8),  t_data(:,8), 100, '*')
             yerr = abs(e_data(:,8)-t_data(:,8));
             errorbar(timedata(:,8),  t_data(:,8), yerr, '*');
             hold on            
             % --- shade 50% CI
             % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',purp,purp,0,0.2); 
             % --- shade 90% CI
             [ph,msg] = jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',purp,purp,1,0.2);
             hold on
             legend('model', 'data', '', 'CI-90%', 'Location', 'NorthWest')
             xlabel('Time (hour)');
             ylabel(speciesNames([30])); 
             ylim([0,1.5]);
             set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Chen et al. (2021)'); t.FontSize = 14;
             hold on
             lgd = legend('model', 'data',  '', 'CI-90%', 'Location', 'SouthEast');
             lgd.FontSize = 14;
             grid on
             
         
          end
     end
    end
    if GLU==0 && LPS==1
        Y_var = Yout(:,[23,24,25,27,28,20]);             % Training variables
    num_outputvar = size(Y_var,2);                   % number of variables
    speciesNames = params{4};                         
    figure(7)
     for i = 1:num_outputvar
         for t=1:length(Time)
             p = kde(Y_var(t,i), 0.05);
             meanholder(t) = meanKDE(p);             % returns mean of KDE
             sdholder(t) = sqrt(covar(p));           % returns variance of KDE
             q = quantile(sample(p,10000),[0.05 0.25 0.75 0.95]); % when sample() is a vector, q(i) is the p(i)th quartiles (0.05 or 0.25 or ...), q is same size as columns of sample()
             qholder(t,i,:) = q;
         end
         % convert any possible negative quantile values to a very small number, 1e-100
         for k = 1:4                                 % determined by the confidences you want to calculate (1:2 if you only want one CI say 90th with quantile 0.05 and 0.95)
             for iS = 1:length(Time)
                 if qholder(iS,i,k) < 0
                     qholderPositive(iS,i,k) = 1e-100;
                 else
                     qholderPositive(iS,i,k) = qholder(iS,i,k);
                 end
             end
         end
         
         if i == 1
             subplot(2,3,i)
             % --- Prediction
             plot(Tout, Y_var(:,i), 'LineWidth', 2, 'color', purp)
             hold on
             % --- Validation data
             scatter(timev(:,1), v_data(:,1), 100, '*')
             yerr = abs(ev_data(:,1) - v_data(:,1));
             errorbar(timev(:,1), v_data(:,1), yerr, '*');
             hold on            
             % --- shade 50% CI
             % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',purp,purp, 0,0.2); 
             hold on
             %--- shade 90% CI
             [ph,msg]=jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',purp,purp,1,0.2);
             hold on            
             xlabel('Time (hour)');
             ylabel(speciesNames([23]));
             ylim([0,1.5]);
             set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Ayala et al. (2019)'); t.FontSize = 14;
             grid on
             
         elseif i == 2             
             subplot(2,3,i)
             % --- Prediction
             plot(Tout, Y_var(:,i), 'LineWidth',2, 'color', purp)
             hold on
             % --- Validation data
             scatter(timev(:,3), v_data(:,3),  100, '*')
             yerr = abs(ev_data(:,3) - v_data(:,3));
             errorbar(timev(:,3), v_data(:,3),  yerr, '*');
             hold on             
             % --- shade 50% CI
             % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',purp,purp,0,0.2); 
             hold on
             % --- shade 90% CI
             [ph,msg]=jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',purp,purp,1,0.2);
             hold on
             xlabel('Time (hour)');
             ylabel(speciesNames([24]));
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);    
             hold on
             t = title('Ayala et al. (2019)'); t.FontSize = 14;         
             grid on
          
         elseif i==3 
             subplot(2,3,i);
             % --- Prediction
             plot(Tout, Y_var(:,i), 'LineWidth', 2, 'color', purp)
             hold on
             % --- Validation data
             scatter(timev(:,2), v_data(:,2),  100, '*')
             yerr = abs(ev_data(:,2)-v_data(:,2));
             errorbar(timev(:,2), v_data(:,2), yerr, '*');
             hold on             
             % --- shade 50% CI
             % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',purp,purp,0,0.2); 
             hold on
             % --- shade 90% CI
             [ph,msg]=jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',purp,purp,1,0.2);
             hold on
             xlabel('Time (hour)');
             ylabel(speciesNames([25])); 
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('He et al. (2020)'); t.FontSize = 14;
             grid on

         elseif i==4
             subplot(2,3,i);
             % --- Prediction
             plot(Tout, Y_var(:,i), 'LineWidth', 2, 'color', purp)
             hold on
             % --- Validation data
             scatter(timev(:,4), v_data(:,4), 100, '*')
             yerr = abs(ev_data(:,4)-v_data(:,4));
             errorbar(timev(:,4), v_data(:,4), yerr,'*');
             hold on
             % --- shade 50% CI
             % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',purp,purp,0,0.2);
             % --- shade 90% CI
             [ph,msg]=jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',purp,purp,1,0.2);
             hold on
             xlabel('Time (hour)');
             ylabel(speciesNames([27])); 
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Yang et al. (2008)'); t.FontSize = 14;
             grid on

         elseif i==5
              subplot(2,3,i)
              % --- Prediction
              plot(Tout, Y_var(:,i), 'LineWidth', 2, 'color', purp);
              hold on
              % --- Validation data
              scatter(timev(:,8), v_data(:,8), 100, '*')
              yerr = abs(ev_data(:,8)-v_data(:,8));
              errorbar(timev(:,8), v_data(:,8), yerr, '*');
              % --- shade 50% CI
              % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',purp,purp,0,0.2); 
              % --- shade 90% CI
              [ph,msg] = jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',purp,purp,1,0.2);
              hold on              
              xlabel('Time (hour)');
              ylabel(speciesNames([28]));
              ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18); 
              hold on
              t = title('Du et al. (2015)'); t.FontSize = 14; 
              hold on
            
              grid on

         elseif i==6
             
              subplot(2,3,i)
              % --- Prediction
              plot(Tout,Y_var(:,i), 'LineWidth', 2, 'color', purp);
              hold on
              % --- Validation data
              scatter(timev(:,7),  v_data(:,7), 100, '*')
              yerr = abs(ev_data(:,7)-v_data(:,7));
              errorbar(timev(:,7),  v_data(:,7), yerr, '*');
              hold on              
              % --- shade 50% CI
              % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',purp,purp,0,0.2); 
              hold on
              % --- shade 90% CI
              [ph,msg]=jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',purp,purp,1,0.2);
              hold on
              xlabel('Time (hour)');
              ylabel(speciesNames([20])); 
              ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
              hold on
              t = title('Ramanathan et al. (2002)'); t.FontSize = 14;
              hold on
              lgd = legend('model', 'data', '', 'CI-90%', 'Location', 'SouthEast');
              lgd.FontSize = 14;
              grid on
         end
          
     end
    end
    if GLU==1 && LPS==1
        Y_var = Yout(:,[23,24,27]);             % Training variables
    num_outputvar = size(Y_var,2);                   % number of variables
    speciesNames = params{4};                         
    figure(8)
     for i = 1:num_outputvar
         for t=1:length(Time)
             p = kde(Y_var(t,i), 0.05);
             meanholder(t) = meanKDE(p);             % returns mean of KDE
             sdholder(t) = sqrt(covar(p));           % returns variance of KDE
             q = quantile(sample(p,10000),[0.05 0.25 0.75 0.95]); % when sample() is a vector, q(i) is the p(i)th quartiles (0.05 or 0.25 or ...), q is same size as columns of sample()
             qholder(t,i,:) = q;
         end
         % convert any possible negative quantile values to a very small number, 1e-100
         for k = 1:4                                 % determined by the confidences you want to calculate (1:2 if you only want one CI say 90th with quantile 0.05 and 0.95)
             for iS = 1:length(Time)
                 if qholder(iS,i,k) < 0
                     qholderPositive(iS,i,k) = 1e-100;
                 else
                     qholderPositive(iS,i,k) = qholder(iS,i,k);
                 end
             end
         end
         
         if i == 1
             subplot(1,3,i)
             % --- Prediction
             plot(Tout, Y_var(:,i), 'LineWidth', 2, 'color', purp)
             hold on
             % --- Validation data
             scatter(timev(:,1), v_data(:,1), 100, '*')
             yerr = abs(ev_data(:,1) - v_data(:,1));
             errorbar(timev(:,1), v_data(:,1), yerr, '*');
             hold on            
             % --- shade 50% CI
             % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',purp,purp, 0,0.2); 
             hold on
             %--- shade 90% CI
             [ph,msg]=jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',purp,purp,1,0.2);
             hold on            
             xlabel('Time (hour)');
             ylabel(speciesNames([23]));
             ylim([0,1.5]);
             set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Ayala et al. (2019)'); t.FontSize = 14;
             grid on
             
         elseif i == 2             
             subplot(1,3,i)
             % --- Prediction
             plot(Tout, Y_var(:,i), 'LineWidth',2, 'color', purp)
             hold on
             % --- Validation data
             scatter(timev(:,3), v_data(:,3),  100, '*')
             yerr = abs(ev_data(:,3) - v_data(:,3));
             errorbar(timev(:,3), v_data(:,3),  yerr, '*');
             hold on             
             % --- shade 50% CI
             % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',purp,purp,0,0.2); 
             hold on
             % --- shade 90% CI
             [ph,msg]=jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',purp,purp,1,0.2);
             hold on
             xlabel('Time (hour)');
             ylabel(speciesNames([24]));
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);   
             hold on
             t = title('Ayala et al. (2019)'); t.FontSize = 14;          
             grid on
          
         else 
             subplot(1,3,i);
             % --- Prediction
             plot(Tout, Y_var(:,i), 'LineWidth', 2, 'color', purp)
             hold on
             % --- Validation data
             scatter(timev(:,4), v_data(:,4), 100, '*')
             yerr = abs(ev_data(:,4)-v_data(:,4));
             errorbar(timev(:,4), v_data(:,4), yerr,'*');
             hold on
             % --- shade 50% CI
             % [ph,msg]=jbfill(Time,qholderPositive(:,i,2)',qholderPositive(:,i,3)',purp,purp,0,0.2);
             % --- shade 90% CI
             [ph,msg]=jbfill(Time,qholderPositive(:,i,1)',qholderPositive(:,i,4)',purp,purp,1,0.2);
             hold on
             xlabel('Time (hour)');
             ylabel(speciesNames([27])); 
             ylim([0,1.5]); set(gca,'FontName','Arial','FontSize',18);
             hold on
             t = title('Yang et al. (2008)'); t.FontSize = 14;
             hold on
             lgd = legend('model', 'data', '', 'CI-90%', 'Location', 'SouthEast');
             lgd.FontSize = 14;
             grid on

         end
          
     end
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
        plot(Tout, Yout(:,var(j)), linest(j), 'LineWidth', 2)
        hold on
        ylabel("Activity"); xlabel('Time (hour)')
        set(gca,'FontName','Arial','FontSize',18);
        hold on
     lgd = legend(params{4}([3,8,10]), 'Location', 'SouthEast');
     lgd.FontSize = 14;

    end
    var = [9,14,15,19];
    for j = 1:length(var)
        subplot(2,3,2)
        plot(Tout, Yout(:,var(j)), linest(j), 'LineWidth', 2); 
        hold on
        ylabel("Activity"); xlabel('Time (hour)')
        set(gca,'FontName','Arial','FontSize',18);
        lgd = legend(params{4}([9,14,15,19]), 'Location', 'SouthEast');
        lgd.FontSize = 14;
        
    end
    var = [7,11,18];
    for j=1:length(var)
        subplot(2,3,3)
        plot(Tout, Yout(:,var(j)), linest(j), 'LineWidth', 2); 
        hold on
        ylabel("Activity"); xlabel('Time (hour)')
        set(gca,'FontName','Arial','FontSize',18);
        lgd = legend(params{4}([7,11,18]), 'Location', 'SouthEast');
        lgd.FontSize = 14;
    end
    var = [4,5,26,16,17];
    for j=1:length(var)
        subplot(2,3,4)
        plot(Tout, Yout(:,var(j)), linest(j), 'LineWidth', 2); 
        hold on
        ylabel("Activity"); xlabel('Time (hour)')
        set(gca,'FontName','Arial','FontSize',18);
        lgd = legend(params{4}([4,5,26,16,17]), 'Location', 'SouthEast');
        lgd.FontSize = 14;
    end
    var = [28,29,30];
    for j=1:length(var)
        subplot(2,3,5)
        plot(Tout, Yout(:,var(j)), linest(j), 'LineWidth', 2)
        hold on
        ylabel("Activity"); xlabel('Time (hour)')
        set(gca,'FontName','Arial','FontSize',18);
        lgd = legend(params{4}([28,29,30]), 'Location', 'SouthEast');
        lgd.FontSize = 14;
    end

end
