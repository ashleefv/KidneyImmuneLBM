function [Tout, Ymean] = coupledODE_IVV_multirun(tspan, y0, params, p_params, mode, state, global_p_best, p_fitted, error_fitted)


global number_ctrl time_ctrl density diameter GC_conc GC_time GC_LB GC_UB time_lee glu_UB glu_LB time_finch glu_finch LB_lee UB_lee glucose_lee

GLU = params{1}(1,1);

indexforNumber = 36;
indexforDiameter = 37;
intv = "none";

blue = 	[0 0.4470 0.7410];
% Time (long-term mice sim.)
    start_time = 2; %weeks
    start_time_h = start_time*7*24;
    end_time = 20; %weeks
    end_time_h = end_time*7*24;
    tspan = start_time_h:1:end_time_h; % hours

if mode == 3
    MC = load('data/MC_25_fen_multirun.mat');
    credible = MC.credible;
    opts = odeset(RelTol=1e-3);
    Nn = 100;
    if state == 'diab_mice'
        glu_sampled = zeros(11,Nn);
        rng("twister") % Default random number generator algorithm with seed = 0 to ensure that we generate the same sequence of draws
        glucose_data_Lee_sd = abs(glucose_lee - LB_lee);
        glucose_data_Finch_sd = abs(glu_finch - glu_UB); 
        subsetIdx = find(time_lee>=6*24*7);
        sigma_data = [glucose_data_Lee_sd(subsetIdx)'  glucose_data_Finch_sd];
%         for Nstep = 1:Nn
%             for i = 1:length(GC_time)
% %                glu_sampled(i,Nstep) = unifrnd(GC_LB(:,i), GC_UB(:,i)); %
%                 glu_sampled(i,Nstep) = normrnd(GC_conc(:,i), sigma_data(i)); %
%             end
        for i = 1:length(GC_time)
%                glu_sampled(i,Nstep) = unifrnd(GC_LB(:,i), GC_UB(:,i)); %
            glu_sampled(i,:) = normrnd(GC_conc(i), sigma_data(i),[1,Nn]); %
        end
        for Nstep = 1:Nn

            [t, y] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,params,p_params, state, glu_sampled(:,Nstep), intv);
            YstepP(Nstep,:,:) = real(y);
            fprintf('run %i finished\n', Nstep)

        end

        Ymean(1,:,:) = mean(YstepP(1:Nn,:,:));
        Ystd(1,:,:) = std(YstepP(1:Nn,:,:));
        Tout = t;
        
        % Glucose plots
        %GLU_p(:,[1:Nn]) = zeros([length(Tout),Nn);
        %GLU_p(1,:) = mean(GC_conc);


        for Nstep=1:Nn
            %Gp = step_function(glu_sampled(:,Nstep));
            for Tt = start_time_h:1:tspan(end)
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

        Ymean0(1,:,:) = mean(YstepP0(1:Nn,:,:));
        Tout = t;
    end

time_g = start_time_h:length(GLU_p(:,1));
   
var = 1:37;

%%

fig = figure(53);
figname = 'FigC';

for i=var
     subplot(5,8,i)
     plot(time_g./(24*7), YstepP(1:Nn,:,i),'LineWidth', 1.2);
    
     hold on
     % for aesthetics make the mean of samples black
     MEANYstepP(:,i)=mean(YstepP(1:Nn,:,i),1);
     plot(time_g./(24*7), MEANYstepP(:,i),'LineWidth', 1.2,'Color','k');
     ylabel(params{4}(i))
     %xlabel('Time (week)')
    set(gca,'FontName','Arial','FontSize',8)
end

% Common x-axis label
han = axes(fig, 'visible', 'off'); 
han.XLabel.Visible = 'on';
xlabel(han, 'Time (weeks)','FontName','Arial','FontSize',8);

widthInches = 9;
heightInches = 5;
run('ScriptForExportingImages.m')   

% reduce YstepP to only the times that correspond to GC_time
[~, idx] = ismember(GC_time, time_g);
valid_idx = idx(idx > 0);
YstepP_GC = YstepP(:, valid_idx, :);
varsofinterest = [1, indexforNumber, indexforDiameter];
% varsofinterest vs. time histograms over Nn samples
for i = varsofinterest
    fignumber = str2num(['530',num2str(i)]);
    fig = figure(fignumber);
    for j=1:length(GC_time)
        subplot(4,3,j)
        histogram(YstepP_GC(:,j,i),'LineWidth', 1.2);
        xlabel(GC_time(j)/24/7)
        set(gca,'FontName','Arial','FontSize',8)
        title(params{4}(i))
    end
    
end

fig = figure(54);
figname = 'FigD';
tiledlayout(length(GC_time),length(varsofinterest))
for j=1:length(GC_time) % rows
    for i = varsofinterest % columns
        nexttile
        histogram(YstepP_GC(:,j,i),'LineWidth', 1.2);
        if i == 1 % common ylabels on the side of each row
            ylabel(['Week ' num2str(GC_time(j)/24/7)],'FontWeight', 'bold')
        end
        if j == 1 % common titles on the top of each column
            title(params{4}(i))
        end
        set(gca,'FontName','Arial','FontSize',8)
    end

end
% for i = varsofinterest
%     mean(YstepP_GC(:,end,i))
%     std(YstepP_GC(:,end,i))
% end

widthInches = 6;
heightInches = 9;
run('ScriptForExportingImages.m')   
% Figure 4 for multiple glucose sample runs
%FINCH figure 2 data for time series
    pop_diameter = load("data\dbmice_diameter_population.csv");
    pop_density = load("data\dbmice_density_population.csv");
    time_g = start_time_h:1:end_time_h;

    % use the FINCH data to determine the standard deviations rather than
    % digitizing the error bars
    pop_time_vector = [6,10,15,20]; % weeks
    pop_diameter_mean = zeros(size(pop_time_vector));
    pop_diameter_sd = zeros(size(pop_time_vector));
    pop_density_mean = zeros(size(pop_time_vector));
    pop_density_sd = zeros(size(pop_time_vector));
    for i = 1:length(pop_time_vector)
        time = pop_time_vector(i);
        matching_rows = pop_diameter(:,1) == time;
        diameters_at_time = pop_diameter(matching_rows, 2);      
        pop_diameter_mean(i) = mean(diameters_at_time);
        pop_diameter_sd(i) = std(diameters_at_time);
        matching_rows = pop_density(:,1) == time;
        density_at_time = pop_density(matching_rows, 2);
        pop_density_mean(i) = mean(density_at_time);
        pop_density_sd(i) = std(density_at_time);
    end

    figure(3) 
    figname = 'Fig3';
    % this figure runs at the output for the mean of the glu_sampled
    % glucose variations
    subplot(2,2,1); box;
    
    hold on 
    % means of control data in FINCH fig. 2
    scatter(pop_time_vector(1), number_ctrl(1), 100, 's', 'filled', 'r');
    hold on 
    % scatter([6,10,15,20], [5.65, 4.50, 4.08, 4.14], 50, 'filled', 'k')
    % hold on
    %errorbar(pop_time_vector, density, abs([5.65, 4.50, 4.08, 4.14] - [6.6, 4.96,4.98,5.41]), abs([5.65, 4.50, 4.08, 4.14] - [4.98, 4.11, 3.15, 2.95]), 'ko', 'MarkerFaceColor', 'k',  'LineWidth', 0.75);
    errorbar(pop_time_vector, density, pop_density_sd, 'ko', 'MarkerFaceColor', 'k',  'LineWidth', 0.75);

    hold on
    scatter(pop_density(:,1), pop_density(:,2), 'o', 'k')
    hold on
    plot(Tout/(24*7), Ymean(1,:,indexforNumber), 'LineWidth', 1.2, 'Color', 'k'); 
    % the standard deviations are within the ODE solver tolerance
    % plot(Tout/(24*7), Ymean(1,:,indexforNumber)+Ystd(1,:,indexforNumber), 'LineWidth', 1.2, 'Color', 'r'); 
    % plot(Tout/(24*7), Ymean(1,:,indexforNumber)-Ystd(1,:,indexforNumber), 'LineWidth', 1.2, 'Color', 'r'); 

    hold on
    [ph,msg] = jbfill(time_g/(24*7), credible(:,1,indexforNumber+1)', credible(:,2,indexforNumber+1)', blue, blue, 1, 0.2);
    xlabel('Time (weeks)');
    ylabel('Fenestration Number');
    xlim([0 21])
    ylim([2 7])
    ax = gca;  ax.FontSize = 8;    
%     t = title(['Finch et al., JASN (2022)']); t.FontSize = 8;

    subplot(2,2,2); box;
    hold on
    
    scatter(pop_time_vector(1), 47.91, 100, 's', 'filled', 'r');
    hold on 
    % scatter([6,10,15,20], [50.74, 60.19, 73.65, 74.63], 50, 'filled', 'k')
    % hold on
    %errorbar(pop_time_vector, diameter', abs([50.74, 60.19, 73.65, 74.63] - [55.12,63.9,80.48,80]), 'ko', 'MarkerFaceColor', 'k', 'LineWidth', 0.75);
    errorbar(pop_time_vector, diameter', pop_diameter_sd, 'ko', 'MarkerFaceColor', 'k', 'LineWidth', 0.75);

    hold on
    scatter(pop_diameter(:,1), pop_diameter(:,2), 'o',  'k')
    hold on
    plot(Tout/(24*7), Ymean(1,:,indexforDiameter), 'LineWidth', 1.2, 'Color', 'k'); 
    % the standard deviations are within the ODE solver tolerance
    % plot(Tout/(24*7), Ymean(1,:,indexforDiameter)+Ystd(1,:,indexforDiameter), 'LineWidth', 1.2, 'Color', 'r'); 
    % plot(Tout/(24*7), Ymean(1,:,indexforDiameter)-Ystd(1,:,indexforDiameter), 'LineWidth', 1.2, 'Color', 'r'); 
    hold on
    [ph,msg] = jbfill(time_g/(24*7), credible(:,1,indexforDiameter+1)', credible(:,2,indexforDiameter+1)', blue, blue, 1, 0.2);
    ylabel('Fenestration Diameter (nm)'); xlabel('Time (weeks)')
%    legend('Model', '95% credible interval', 'Control data', 'Diabetes data', 'Diabetes population')
    ax = gca;  ax.FontSize = 8;    
    xlim([0 21])
    ylim([40 100])
%     t = title(['Finch et al., JASN (2022)']); t.FontSize = 8;

    labelstring = {'A', 'B'};
    for v = 1:2
        subplot(2,2,v)
        hold on
        text(-0.15, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize',8)
        set(gca,'FontName','Arial','FontSize',8)
    end

        % overall legend
    legend( 'Control data', 'Diabetes data (mean+/-SD)',  'Diabetes data (individual)','Model', '95% credible interval')
    h = legend('Location','southoutside', 'Orientation', 'horizontal');
    h.NumColumns = 3;
    p = [0.5 0.45 0.03 0.03];
    set(h,'Position', p,'Units', 'normalized');

widthInches = 6;
heightInches = 4.6;
run('ScriptForExportingImages.m') 


end
%% Parameter Variability

if mode == 4

global_p_best = csvread("recal-param\fen_combined_25.csv");
fitted_p = csvread("recal-param\fen_combined_fitted_25.csv");
p_fitted = fitted_p(:,1:end-1);
error_fitted = fitted_p(:,end);


GB = params;
pPh = p_params;


pPh(4) = global_p_best(:,1); % ks
pPh(5) = global_p_best(:,2); % kd
GB{2}(33) = global_p_best(:,3); % tau 33
pPh(6) = global_p_best(:,4); % kloss
pPh(7) = global_p_best(:,5); % N_ss2

% sort parameters by minimum error
Params_error = [error_fitted p_fitted]; % coef_fitted is the parameter output from all the LHS fmincon runs
error_column = 1; % first column stores error_fitted
[sorted_Params_error, index] = sortrows(Params_error,error_column);
param_posterior = sorted_Params_error(:,2:end);


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

p_posterior = sorted_Params_error(accep_id,2:end); % acceptable parameters within 1.2*min(SSE)
Ns = 100; % number of samples

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
    dp{2}(33) = population(j,3);
    pPh(6) = population(j,4);
    pPh(7) = population(j,5);
    
    options = [];
    [Tout, Yout] = coupledODE_IVV_run(tspan, y0, dp, pPh, 0, state, GC_conc');
    Y_param_var(j,:,:) = Yout(:,:);
    

end



for j = 1:length(params{2}(:))
    for t = 1:length(Tout)
        credible(t,:,j) = quantile(abs(Y_param_var(:,t,j)), [0.025 0.975]);  % diameter
        Yp_mean(t,j) = mean(abs(Y_param_var(:,t,j)));
    end
%save('data/MC_25_fen_multirun.mat')
%note that this was last saved when we had one more variable in a legacy
%model. Rather than rerunning this time-intensive MC simulation, we simple
%reuse the previously generated credible intervals matrix credible from 
%data/MC_25_fen_multirun.mat and offset indexforDiameter and indexforNumber 
%to account for the legacy variable, which doesn't impact the calculations 
%for anything else in the model. 
end

%%
[Tbest, Ybest] = coupledODE_IVV_run(tspan, y0, GB, pPh, 0, state, GC_conc'); % simulated at best fit parameters

% individual mice data from Finch et al.

pop_diameter = load("data\dbmice_diameter_population.csv");
pop_density = load("data\dbmice_density_population.csv");


figure(101)
time_g = tspan;

hold on
plot(time_g/(24*7), Ybest(:,indexforDiameter), 'color', 'k', 'LineWidth', 1)             % mean of acceptable estimates
hold on
[ph,msg] = jbfill(time_g/(24*7), credible(:,1,indexforDiameter+1)', credible(:,2,indexforDiameter+1)', blue, blue, 1, 0.2);
hold on
scatter(6, 47.91, '*', 'r');
hold on 
scatter([6,10,15,20], [50.74, 60.19, 73.65, 74.63], 50, 'filled', 'k')
hold on
errorbar([6,10,15,20], [50.74, 60.19, 73.65, 74.63], abs([50.74, 60.19, 73.65, 74.63] - [55.12,63.9,80.48,80]), abs([50.74, 60.19, 73.65, 74.63] - [55.12,63.9,80.48,80]), 'o', 'Color',[1  0  0]);
hold on
scatter(pop_diameter(:,1), pop_diameter(:,2), '^', 'k')
ylabel('Fenestration Diameter or Width (nm)'); xlabel('Time (week)')
legend('Model', '95% credible interval', 'control data', 'diabetes data', '', 'diabetes population')

 
figure(102)
plot(time_g/(24*7), Ybest(:,indexforNumber), 'color', 'k', 'LineWidth', 1)             % mean of acceptable estimates
hold on
[ph,msg] = jbfill(time_g/(24*7), credible(:,1,indexforNumber+1)', credible(:,2,indexforNumber+1)', blue, blue, 1, 0.2);
hold on 
scatter(6, 6.3, '*', 'r');
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