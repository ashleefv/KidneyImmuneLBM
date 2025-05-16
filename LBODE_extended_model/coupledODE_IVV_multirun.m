function [Tout, Ymean] = coupledODE_IVV_multirun(tspan, y0, params, p_params, mode, state, global_p_best, p_fitted, error_fitted)


global number_ctrl time_ctrl density diameter GC_conc GC_time GC_LB GC_UB time_lee glu_UB glu_LB time_finch glu_finch LB_lee UB_lee glucose_lee

GLU = params{1}(1,1);
blue = 	[0 0.4470 0.7410];

intv = "none";

if mode == 3
    MC = load('data/MC_25_fen_multirun.mat');
    credible = MC.credible;
    opts=[];
    Nn = 25;
    if state == 'diab_mice'
        
        glu_sampled = zeros(11,Nn);

        for Nstep = 1:Nn
            for i = 1:length(GC_time)
                glu_sampled(i,Nstep) = unifrnd(GC_LB(:,i), GC_UB(:,i)); %

            end

            [t, y] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,params,p_params, state, glu_sampled(:,Nstep), intv);
            YstepP(Nstep,:,:) = real(y);
            fprintf('run %i finished\n', Nstep)

        end

        Ymean(1,:,:) = mean(YstepP([1:Nn],:,:));
        Tout = t;
        
        % Glucose plots
        %GLU_p(:,[1:Nn]) = zeros([length(Tout),Nn);
        %GLU_p(1,:) = mean(GC_conc);


        for Nstep=1:Nn
            %Gp = step_function(glu_sampled(:,Nstep));
            for Tt = [336:1:tspan(end)]
                if (Tt >= time_lee(1) && Tt <= time_lee(5))
                    GLU_p(Tt,Nstep) =  0.051*(Tt)  - 9.38;
                else
                    GLU_p(Tt,Nstep) = step_function(Tt, glu_sampled(:,Nstep)); %double(Gp([Tt]));
                end
            end
        end
    
    elseif state == 'norm_mice'
        for Nstep = 1:Nn
            
            [t, y] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,params,p_params, state, [], intv);
            YstepP0(Nstep,:,:) = real(y);
            fprintf('run %i finished\n', Nstep)

        end

        Ymean0(1,:,:) = mean(YstepP0([1:Nn],:,:));
        Tout = t;
    end

time_g = [336:length(GLU_p(:,1))];
%%

% A = 1; %2.06; %p_params(end);
% B = 1;
% 
% 
% FenC(1) =  y0(37); 
% 
% for Tt = [2:1:length(tspan)]
%        
% 
%        FenC(Tt) = ((1/(1 + A*Ymean(1,Tt,2))) + (B*Ymean(1,Tt,31)/(1 + B*Ymean(1,Tt,31))) - 2*(1/(1 + A*Ymean(1,Tt,2)))*(B*Ymean(1,Tt,31)/(1 + B*Ymean(1,Tt,31)))); 
%        %wt = exp((Tt)/2500);
%        
% end
% Ymean(1,:,37) = FenC(:); 

% pred_in = interp1(Tout, Ymean(1,:,37), invivo_time([12,16,32,36],11)*24);
% error(1,1) = sum((pred_in - invivo_train([12,16,32,36],11)).^2);
%%

% 
% 
% figure(2)
% subplot(1,2,2)
% hold on
% plot(time_g/(24*7), Ymean(1,:,1), 'LineWidth', 3, 'Color', 'k');
% xlabel('Time (weeks)'); %xlim([0,22]);
% ylabel('Normalized Glucose'); %ylim([0,40])
% ax = gca; ax.FontSize = 20;
% hold on
% ax = gca; ax.FontSize = 20;
% 
% figure(2)
% subplot(1,2,1)
% plot(time_g/(24*7), mean(GLU_p(time_g, [1:Nn]),2), 'LineWidth', 3, 'Color', 'k');
% hold on; scatter(time_finch/(7*24), glu_finch, 'MarkerFaceColor',[1 0 0],'Marker','^','Color',[1 0 0])
% hold on; errorbar(time_finch/(7*24), glu_finch, abs(glu_finch - glu_LB), abs(glu_finch - glu_UB), '^', 'Color',[1  0  0]);
% hold on; scatter(time_lee/(7*24), glucose_lee, 'MarkerFaceColor',[0  0  1],'Marker','o','Color',[0  0  1])
% hold on; errorbar(time_lee/(7*24), glucose_lee, abs(glucose_lee - LB_lee), abs(glucose_lee - UB_lee), 'o', 'Color',[0  0  1]);
% xlabel('Time (weeks)'); xlim([0,21]);
% ylabel('Glucose (mmol/l)'); ylim([0,55]);
% 
% legend('GLU_p, sampled', 'Finch et al. (2022)', '', 'Lee et al. (2018)', '');
% ylim([0,50]);
% ax = gca; ax.FontSize = 20;



 
    pop_diameter = load("data\dbmice_diameter_population.csv");
    pop_density = load("data\dbmice_density_population.csv");
    time_g = [336:1:3360];

    
    figure(3) 
    subplot(1,2,1); box;
    hold on

    plot(Tout/(24*7), Ymean(1,:,[38]), 'LineWidth', 1.2, 'Color', 'k'); 
    hold on
    hold on
    [ph,msg] = jbfill(time_g/(24*7), credible(:,1,38)', credible(:,2,38)', blue, blue, 1, 0.2);
    hold on
    scatter([6], [47.91], 100, 's', 'filled', 'r');
    hold on 
    scatter([6,10,15,20], [50.74, 60.19, 73.65, 74.63], 50, 'filled', 'k')
    hold on
    errorbar([6,10,15,20], [50.74, 60.19, 73.65, 74.63], abs([50.74, 60.19, 73.65, 74.63] - [55.12,63.9,80.48,80]), abs([50.74, 60.19, 73.65, 74.63] - [55.12,63.9,80.48,80]), 'o', 'Color','k', 'LineWidth', 0.75);
    hold on
    scatter(pop_diameter(:,1), pop_diameter(:,2), 'o',  'k')
    ylabel('Fenestration Diameter (nm)'); xlabel('Time (week)')
%    legend('Model', '95% credible interval', 'Control data', 'Diabetes data', '', 'Diabetes population')
    ax = gca;  ax.FontSize = 18;    
%     t = title(['Finch et al., JASN (2022)']); t.FontSize = 14;



    figure(3) 
    subplot(1,2,2); box;
    hold on
    plot(Tout/(24*7), Ymean(1,:,37), 'LineWidth', 1.2, 'Color', 'k'); 
    [ph,msg] = jbfill(time_g/(24*7), credible(:,1,37)', credible(:,2,37)', blue, blue, 1, 0.2);
    hold on 
    scatter([6], [6.3], 100, 's', 'filled', 'r');
    hold on 
    scatter([6,10,15,20], [5.65, 4.50, 4.08, 4.14], 50, 'filled', 'k')
    hold on
    errorbar([6,10,15,20], [5.65, 4.50, 4.08, 4.14], abs([5.65, 4.50, 4.08, 4.14] - [6.6, 4.96,4.98,5.41]), abs([5.65, 4.50, 4.08, 4.14] - [4.98, 4.11, 3.15, 2.95]), 'o', 'Color', 'k', 'LineWidth', 0.75);
    hold on
    scatter(pop_density(:,1), pop_density(:,2), 'o', 'k')
    
    xlabel('Time (week)');
    ylabel('Fenestration Number');
    legend('Model', '95% credible interval', 'Control data', 'Diabetes data (mean)', '', 'Diabetes data (individual)')
    ax = gca;  ax.FontSize = 18;    
%     t = title(['Finch et al., JASN (2022)']); t.FontSize = 14;

    hold off;
         
%%
if state == 'diab_mice'
    
    figure(20); 
    meanGLP = mean(GLU_p([336:3360],:),2);
    yyaxis left
    % plot([336:length(GLU_p(:,:))]/(24*7), GLU_p([336:3360],:), 'LineWidth', 1, 'Color', [.7 .7 .7]); 
    % jbfill([336:3360]/24/7, min(GLU_p([336:3360],:)'), max(GLU_p([336:3360],:)'), [.7 .7 .7], [.7 .7 .7], 1, 0.2)
    plot(Tout/(24*7), meanGLP, 'LineWidth', 1, 'LineStyle' , '--', 'Color', [0 0.4470 0.7410]); 
    hold on; scatter(time_finch/(7*24), glu_finch, 'MarkerFaceColor',[1 0 0],'Marker','^','Color',[1 0 0])
    hold on; errorbar(time_finch/(7*24), glu_finch, abs(glu_finch - glu_LB), abs(glu_finch - glu_UB), '^', 'Color',[1  0  0]);
    hold on; scatter(time_lee/(7*24), glucose_lee, 'MarkerFaceColor',[0  0  1],'Marker','o','Color',[0  0  1])
    hold on; errorbar(time_lee/(7*24), glucose_lee, abs(glucose_lee - LB_lee), abs(glucose_lee - UB_lee), 'o', 'Color',[0  0  1]);
    xlabel('Time (weeks)'); xlim([0,21]);
    ylabel('Glucose (mmol/l)'); 
    ax = gca; ax.FontSize = 20;
    hold on
    legend('Model', 'Finch et al. (2022)', '', 'Lee et al. (2018)', '', 'Normalized Glucose');
    ylim([7.07,52]);
    x0=10;
    y0=10;
    width=800;
    height=800;
    set(gcf,'position',[x0,y0,width,height])
    ax = gca; ax.FontSize = 20;

    
    % norm_glu(:) = (glu_finch(:) - 10) / (max(glu_UB) - 10);
    % norm_UB(:) = (glu_UB(:) - 10)/(max(glu_UB) - 10);
    % norm_LB(:) = (glu_LB(:) - 10)/(max(glu_UB) - 10);

    yyaxis right
    plot(Tout/(24*7), Ymean(1,:,1), 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
    
    ylabel('Normalized Glucose');% xlim([0,22]); 
    box on

    xlabel('Time (weeks)');
    ylabel('Normalized Glucose'); ylim([0,1.6]);
    
else 
    figure(20)
    plot(Tout(time_in)/(24*7), GLU_p(time_in,1), 'LineWidth', 2);
    hold on; scatter(time_finch/(7*24), ctrl_finch, 'MarkerFaceColor',[0 0 0],'Marker','o','Color',[0.75 0.5 0.25])
    hold on; errorbar(time_finch/(7*24), ctrl_finch, abs(ctrl_finch - ctrl_LB), abs(ctrl_finch - ctrl_UB), '*');
    xlabel('Time (weeks)'); %xlim([0,22]);
    ylabel('Glucose (mmol/l)'); %ylim([0,40])
    ax = gca; ax.FontSize = 20;
    hold on
    legend('Model', 'Data', '')
    x0=10;
    y0=10;
    width=800;
    height=800;
    set(gcf,'position',[x0,y0,width,height])
end

% gap width node 30 removed
var = [1:29, 31:36];

% figure(51)
% for i=var
%     if i<=29
%      subplot(4,9,i)
%      plot(Tout/(24*7), Ymean(1,:,i),'LineWidth', 1.2, 'Color', 'k');
% 
%      hold on
%     end
%     if i>=31
%      subplot(4,9,i-1)
%      plot(Tout/(24*7), Ymean(1,:,i),'LineWidth', 1.2, 'Color', 'k');
% 
%      hold on
%     end
% 
%      ylabel(params{4}(i))
%      xlabel('Time (week)')
%      xlim([0,21]);
%      ax = gca; ax.FontSize = 14;
% end
% 
% hold off;




%%

figure(52)
for i=var
    if i<=29
     subplot(5,8,i)
     plot(time_g/(24*7), YstepP(:,:,i),'LineWidth', 1.2);
     %M(i) = (Ymean(1,end,i) - Ymean(1,1,i));
    end
    if i>=31
     subplot(5,8,i-1)
     plot(time_g/(24*7), YstepP(:,:,i),'LineWidth', 1.2);
     %M(i) = (Ymean(1,end,i) - Ymean(1,1,i));
    end 
%     
     hold on
     ylabel(params{4}(i))
     xlabel('Time (week)')
end

hold off;


%% Calculate p-value between s.s. protein activity of control and diabetic mice
% YstepP0: simulated output for control mice model
% YstepP: simulated output for diabetic mice model

% state = 'norm_mice';
% for Nstep = 1:Nn
%     
%     [t, y] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,params,p_params, state, [], intv);
%     YstepP0(Nstep,:,:) = real(y);
%     fprintf('run %i finished\n', Nstep)
% 
% end
% 
% r = randsample([1:Nn],Nn-1);
% 
% % HM = mean([YstepP(r,end,[2,31:36]) - YstepP0(r,end,[2,31:36])]));
% % heatmap(HM)
% test_sample(:,:,[1:7]) = [YstepP0(r,end,[2,31:36]), YstepP(r,end,[2,31:36])];
% for i=1:7
%     [p, tbl, stats] = anova1(test_sample(:,:,i));
% end


%%%[hypo,pvalue,ci,stats] = ttest(YstepP(:,end,33), 0, 'Alpha', 0.05);
end
%% Parameter Variability

if mode == 4

global_p_best = csvread("recal-param\fen_combined_25.csv");
fitted_p = csvread("recal-param\fen_combined_fitted_25.csv");
p_fitted = fitted_p(:,[1:end-1]);
error_fitted = fitted_p(:,end);

% tau_index = [34, 38]; W_index = [43, 47]; n_index = []; k_index = [43, 47]; 
% size_tau = size(tau_index,2);
% size_k = size(k_index,2); % EC50
% size_W = size(W_index, 2);

GB = params;
pPh = p_params;
% GB{2}(tau_index) = global_p_best(1:size_tau);
% GB{1}(1,W_index) = global_p_best(size_tau+1:size_W+size_tau);
% GB{1}(3,k_index) = global_p_best(size_tau+size_W+1:size_W+size_k+size_tau);


pPh(4) = global_p_best(:,1); % kd
pPh(5) = global_p_best(:,2); % ke
GB{2}(34) = global_p_best(:,3); % tau 34
pPh(6) = global_p_best(:,4); % kloss
pPh(7) = global_p_best(:,5); % N_ss2

% pPh(1) = global_p_best(:,1); % N_ss1
% pPh(2) = global_p_best(:,2); % nf
% pPh(3) = global_p_best(:,3); % kform


% sort parameters by minimum error
col = [0 0.4470 0.7410];
Params_error = [error_fitted p_fitted]; % coef_fitted is the parameter output from all the LHS fmincon runs
error_column = 1; % first column stores error_fitted
[sorted_Params_error, index] = sortrows(Params_error,error_column);
param_posterior = sorted_Params_error(:,[2:end]);


% Sum of Squared Error (SSE) plot
Params_error = [error_fitted p_fitted]; % coef_fitted is the parameter output from all the LHS fmincon runs
error_column = 1; % first column stores error_fitted
[sorted_Params_error, index] = sortrows(Params_error,error_column);
accep_id = find(sorted_Params_error(:,1) < 1.2*min(sorted_Params_error(:,1)));

figure(100)
%-- generate figure to show the error of all the parameter values from LHS
for i = 1:size(p_fitted, 2)
    plot(1:length(index),sorted_Params_error(:,error_column),'o')
    hold on
    xlabel('number of runs')
    ylabel('sum of squared error (sse)')
end

%%
% Monte Carlo Method
disp('Starting Monte Carlo Simulations ...')

p_posterior = sorted_Params_error(accep_id,[2:end]); % acceptable parameters within 1.2*min(SSE)
Ns = 1000; % number of samples

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


%%
for j = 1:Ns
    dp = params;
    pPh = p_params;
%     dp{2}(tau_index) = population(j,1:size_tau);
%     dp{1}(1,W_index) = population(j,size_tau+1:size_W+size_tau);
%     dp{1}(3,k_index) = population(j,size_tau+size_W+1:size_tau+size_k+size_W);
    pPh(4) = population(j,1);
    pPh(5) = population(j,2);
    dp{2}(34) = population(j,3);
    pPh(6) = population(j,4);
    pPh(7) = population(j,5);
    
    options = [];
    [Tout, Yout] = coupledODE_IVV_run(tspan, y0, dp, pPh, 0, state, GC_conc');
    Y_param_var(j,:,:) = Yout(:,:);
    

end

% yd = Y_param_var(:,:,38);
% revhill = (yd*0.25./(1-yd)).^0.5;
% RelScore_fd = (revhill*47.9 + 47.9);

% %
% Neg = -(Y_param_var(:,:,37)*(0.543) - 0.543);
% RevHill = (Neg.*0.25./(1-Neg)).^0.5;
% ChangeTh = -1*RevHill;
% RelScore = ChangeTh*6.4 + 6.4;





blue = 	[0 0.4470 0.7410];
for j = 1:length(params{2}(:))
    for t = 1:length(Tout)
        credible(t,:,j) = quantile(abs(Y_param_var(:,t,j)), [0.025 0.975]);  % diameter
        %credible(t,:) = quantile(abs(RelScore_fd(:,t)), [0.025 0.975]);  % diameter
        % mean_gap(t,1) = mean(abs(Y_param_var(:,t,30)));
        Yp_mean(t,j) = mean(abs(Y_param_var(:,t,j)));
        %Yp_mean(t,38) = mean(abs(RelScore_fd(:,t)));
    end

end

%%
[Tbest, Ybest] = coupledODE_IVV_run(tspan, y0, GB, pPh, 0, state, GC_conc'); % simulated at best fit parameters

% individual mice data from Finch et al.

pop_diameter = load("data\dbmice_diameter_population.csv");
pop_density = load("data\dbmice_density_population.csv");



figure(101)
time_g = [336:1:3360];

hold on
plot(time_g/(24*7), Ybest(:,38), 'color', 'k', 'LineWidth', 1)             % mean of acceptable estimates
hold on
[ph,msg] = jbfill(time_g/(24*7), credible(:,1,38)', credible(:,2,38)', blue, blue, 1, 0.2);
hold on
scatter([6], [47.91], '*', 'r');
hold on 
scatter([6,10,15,20], [50.74, 60.19, 73.65, 74.63], 50, 'filled', 'k')
hold on
errorbar([6,10,15,20], [50.74, 60.19, 73.65, 74.63], abs([50.74, 60.19, 73.65, 74.63] - [55.12,63.9,80.48,80]), abs([50.74, 60.19, 73.65, 74.63] - [55.12,63.9,80.48,80]), 'o', 'Color',[1  0  0]);
hold on
scatter(pop_diameter(:,1), pop_diameter(:,2), '^', 'k')
ylabel('Fenestration Diameter or Width (nm)'); xlabel('Time (week)')
legend('Model', '95% credible interval', 'control data', 'diabetes data', '', 'diabetes population')


figure(102)
plot(time_g/(24*7), Ybest(:,37), 'color', 'k', 'LineWidth', 1)             % mean of acceptable estimates
hold on
[ph,msg] = jbfill(time_g/(24*7), credible(:,1,37)', credible(:,2,37)', blue, blue, 1, 0.2);
hold on 
scatter([6], [6.3], '*', 'r');
hold on 
scatter([6,10,15,20], [5.65, 4.50, 4.08, 4.14], 50, 'filled', 'k')
hold on
errorbar([6,10,15,20], [5.65, 4.50, 4.08, 4.14], abs([5.65, 4.50, 4.08, 4.14] - [6.6, 4.96,4.98,5.41]), abs([5.65, 4.50, 4.08, 4.14] - [4.98, 4.11, 3.15, 2.95]), 'o');
hold on
scatter(pop_density(:,1), pop_density(:,2), '^', 'k')

xlabel('Time (week)');
ylabel('Fenestration Number');

legend('Model', '95% credible interval', 'control data', 'diabetes data', '', 'diabetes population')
toc
end



%% In vivo  glucose and feedback ON


end