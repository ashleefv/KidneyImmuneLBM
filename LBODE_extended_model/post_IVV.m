function [s_FD_Ym, s_FD_W] = post_IVV(params, y0, tspan, p_params, state, task, Tstop)

global number_ctrl time_ctrl density diameter GC_conc GC_time GC_LB GC_UB time_lee glu_UB glu_LB time_finch glu_finch LB_lee UB_lee glucose_lee

% Time (long-term mice sim.)
    start_time = 2; %weeks
    start_time_h = start_time*7*24;
    end_time = 20; %weeks
    end_time_h = end_time*7*24;
    tspan = start_time_h:1:end_time_h; % hours
indexforNumber = 36;
indexforDiameter = 37;

intv = 'none';

RP = length(params{1}(1,:));
SP = length(params{3}(:));

s_FD_Ym = []; s_FD_W = [];
%% 1: test_knockout
if task == 1
    opts = [];%odeset(AbsTol=1e-8);
    Nn = 100;
    state = 'diab_mice';
    
    glu_sampled = zeros(11,Nn);
    rng("twister") % Default random number generator algorithm with seed = 0 to ensure that we generate the same sequence of draws
    glucose_data_Lee_sd = abs(glucose_lee - LB_lee);
    glucose_data_Finch_sd = abs(glu_finch - glu_UB); 
    subsetIdx = find(time_lee>=6*24*7);
    sigma_data = [glucose_data_Lee_sd(subsetIdx)'  glucose_data_Finch_sd];

%     for Nstep = 1:Nn
%         for i = 1:length(GC_time)
% %             glu_sampled(i,Nstep) = unifrnd(GC_LB(:,i), GC_UB(:,i)); %            
%             glu_sampled(i,Nstep) = normrnd(GC_conc(:,i), sigma_data(i)); %
%         end
    for i = 1:length(GC_time)
%                glu_sampled(i,Nstep) = unifrnd(GC_LB(:,i), GC_UB(:,i)); %
        glu_sampled(i,:) = normrnd(GC_conc(i), sigma_data(i),[1,Nn]); %
    end
    for Nstep = 1:Nn
        [t, y] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,params,p_params, state, glu_sampled(:,Nstep), intv);
        YstepP(Nstep,:) = real(y(end,:));
        %disp(size(YstepP))
    end
    s_FC(:,1,:) = YstepP(:,:); % no treatment case

    %% read data sets for healthy and diabetes at 20 weeks with the samples from Finch
    % FINCH figure 1 E & F
    s_ref=readmatrix('data/FINCH_FENESTRATION_20wk.csv');
    s_ctrl_n = s_ref(:,2); % number data for healthy case
    s_ctrl_d = s_ref(:,3); % diameter data for healthy case
    s_ctrl_db_n = s_ref(1:8,4); % number data for diabetic case
    s_ctrl_db_d = s_ref(1:8,5); % diameter data for diabetic case
    % s_ctrl_n([1:Nn],1,1) = 5.8; % reference is set for initial data value for healthy case
    % s_ctrl_d([1:Nn],1,1) = 50;  % reference is set for initial data value for healthy case
    test_i = [29,32,31,34,2]; % inhibitors: KN93, ML7, Y27632, CalA, CytB
    z_params = params;
    % we plot 3 additional columns before test_i: healthy data, diabetes
    % data, no treatment

    % These are the posterior distributions of the output fenestration number and diameter based on MC sampling of the parameter posteriors    
    MC = load('data/MC_25_fen_multirun.mat');
    Y_param_var = MC.Y_param_var;
    MC_samplesn = Y_param_var(:,end,indexforNumber+1);
    MC_samplesd = Y_param_var(:,end,indexforDiameter+1);

    % run t-tests for comparison of number and diameter to each of these
    % three columns
    testColumns = 3;
    p_c_n = zeros(1,length(test_i)+testColumns); 
    p_c_d = zeros(1,length(test_i)+testColumns);
    p_n_data = zeros(1,length(test_i)+testColumns);
    p_d_data = zeros(1,length(test_i)+testColumns);
    p_n = zeros(1,length(test_i)+testColumns);
    p_d = zeros(1,length(test_i)+testColumns);
    offset = 2;
    % compare diabetes data to control data and store in the offset column
    %p-values for number relative to the healty condition no treatment simulation
    [h, p_c_n(1,offset)] = ttest2(s_ctrl_db_n, s_ctrl_n, 'Alpha', 0.05,'Vartype','unequal');
    %p-values for diameter relative to the healty condition no treatment
    [h, p_c_d(1,offset)] = ttest2(s_ctrl_db_d, s_ctrl_d, 'Alpha', 0.05,'Vartype','unequal');
    
% We have two options for a distribution of samples of simulated untreated
% fenestration number and diameter. One of the MC_samples* were * is n or d
% for number and diameter. The second option selected here is the no
% treatment case sampled over glucose inputs. The latter is closer to the
% normal distribution that is expected in a t-test, so is used here.
refForModeln =  s_FC(:,1,indexforNumber);
refForModeld =  s_FC(:,1,indexforDiameter);
    inh = 0; % no treatment case
        %p-values for number relative to the healty condition no treatment simulation
        [h, p_c_n(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforNumber), s_ctrl_n, 'Alpha', 0.05,'Vartype','unequal');
        %p-values for diameter relative to the healty condition no treatment
        [h, p_c_d(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforDiameter), s_ctrl_d, 'Alpha', 0.05,'Vartype','unequal');

         %p-values for number relative to the diseased condition no treatment DATA
        [h, p_n_data(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforNumber), s_ctrl_db_n, 'Alpha', 0.05,'Vartype','unequal');
        %p-values for diameter relative to the diseased condition no treatment DATA
        [h, p_d_data(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforDiameter), s_ctrl_db_d, 'Alpha', 0.05,'Vartype','unequal');
        [h, p_n(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforNumber), refForModeln, 'Alpha', 0.05,'Vartype','unequal');
        %p-values for diameter relative to the diseased condition no treatment simulation
        [h, p_d(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforDiameter), refForModeld, 'Alpha', 0.05,'Vartype','unequal');
        
       

    for inh = 1:length(test_i)
        % treatment knocks out a pathway via its parameter
        z_params{3}(test_i(inh)) = 0;
        

        for Nstep = 1:Nn
            
            [tn, yn] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,z_params,p_params, state, glu_sampled(:,Nstep), intv);
            YstepP_new(Nstep,:) = real(yn(end,:));
            %fprintf('run %i finished\n', Nstep)

        end
        fprintf('finished treatment %i \n', inh)
       
       
        s_FC(:,inh+1,:) = YstepP_new(:,:);

        %p-values for number relative to the healty condition no treatment simulation
        [h, p_c_n(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforNumber), s_ctrl_n, 'Alpha', 0.05,'Vartype','unequal');
        %p-values for diameter relative to the healty condition no treatment
        [h, p_c_d(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforDiameter), s_ctrl_d, 'Alpha', 0.05,'Vartype','unequal');
         %p-values for number relative to the diseased condition no treatment DATA
        [h, p_n_data(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforNumber), s_ctrl_db_n, 'Alpha', 0.05,'Vartype','unequal');
        %p-values for diameter relative to the diseased condition no treatment DATA
        [h, p_d_data(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforDiameter), s_ctrl_db_d, 'Alpha', 0.05,'Vartype','unequal');
        [h, p_n(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforNumber), refForModeln, 'Alpha', 0.05,'Vartype','unequal');
        %p-values for diameter relative to the diseased condition no treatment simulation
        [h, p_d(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforDiameter), refForModeld, 'Alpha', 0.05,'Vartype','unequal');
        

      % %p-values for number relative to the diseased condition no treatment
      %   [h, p_n(1,inh+1)] = ttest2(s_FC(:,inh+1,indexforNumber), s_FC(:,1,indexforNumber), 'Alpha', 0.05,'Vartype','unequal');
      %   %p-values for diameter relative to the diseased condition no treatment
      %   [h, p_d(1,inh+1)] = ttest2(s_FC(:,inh+1,indexforDiameter), s_FC(:,1,indexforDiameter), 'Alpha', 0.05,'Vartype','unequal');
      %   %p-values for number relative to the healty condition no treatment
      %   [h, p_c_n(1,inh)] = ttest2(s_FC(:,inh+1,indexforNumber), s_ctrl_n, 'Alpha', 0.05,'Vartype','unequal');
      %   %p-values for diameter relative to the healty condition no treatment
      %   [h, p_c_d(1,inh)] = ttest2(s_FC(:,inh+1,indexforDiameter), s_ctrl_d, 'Alpha', 0.05,'Vartype','unequal');


%reset the params before the next knockout
        z_params = params;
        
    end
p_c_n
p_c_d
    p_n_data
p_d_data
    p_n
p_d

    s_FC_d = squeeze(s_FC(:,:,indexforDiameter));
    s_FC_n = squeeze(s_FC(:,:,indexforNumber));

    CI_diameter_ub = mean(s_FC_d,1) + std(s_FC_d,1)*1.96; % upper-bound for diameter predicted
    CI_diameter_lb = mean(s_FC_d,1) - std(s_FC_d,1)*1.96;

    CI_number_ub = mean(s_FC_n,1) + std(s_FC_n,1)*1.96; % upper-bound for number predicted
    CI_number_lb = mean(s_FC_n,1) - std(s_FC_n,1)*1.96; 

%     disp(mean(s_FC_d,1)); disp(std(s_FC_d,1));
%     disp(mean(s_FC_n,1)); disp(std(s_FC_n,1));
% 
% %%
%      s_FC([1:Nn],:,indexforNumber)
%      size(s_FC)
     figure(6)
     figname = 'Fig6';
     subplot(2,1,1); b = bar([mean(s_ctrl_n),mean(s_ctrl_db_n), mean(refForModeln), mean(s_FC(1:Nn,2:end,indexforNumber))], 'white'); xticks(1:length(test_i)+3); 
    xticklabels({'Healthy Data', 'Diabetes Data','No Treatment', 'KN93', 'ML7', 'Y27632', 'CalA', 'CytB'}); 
     ylabel('Fenestration Number');
     b.FaceColor = 'flat';
     b.CData(1,:) = [0 0 1];
     b.CData(2,:) = [0 0 0];
     b.CData(3,:) = [.7 .7 .7];
     for inh = 1:length(test_i)
        b.CData(3+inh,:) = [1 1 1]; 
     end
   hold on;  
   %er = errorbar([2:7], [mean(s_FC([1:Nn],:,indexforNumber))], (CI_number_lb-CI_number_ub)); er.Color = 'r';          er.LineStyle = 'none'; er.LineWidth=2;  
   er = errorbar(1, mean(s_ctrl_n), std(s_ctrl_n)); er.Color = 'r';          er.LineStyle = 'none'; er.LineWidth=2; 
   er = errorbar(2, mean(s_ctrl_db_n), std(s_ctrl_db_n)); er.Color = 'r';          er.LineStyle = 'none'; er.LineWidth=2;
   er = errorbar(3, mean(refForModeln), std(refForModeln)); er.Color = 'r';          er.LineStyle = 'none'; er.LineWidth=2;    
   er = errorbar(4:8, mean(s_FC(1:Nn,2:end,indexforNumber)), std(s_FC(1:Nn,2:end,indexforNumber))); er.Color = 'r';          er.LineStyle = 'none'; er.LineWidth=2;  
   ylim([0 10.5])
   
     callsigstar(2,p_n_data,'k')
     callsigstar(1,p_c_n,'b')

     set(gca,'FontSize',8)

     figure(6); subplot(2,1,2); B = bar([mean(s_ctrl_d), mean(s_ctrl_db_d), mean(refForModeld), mean(s_FC(1:Nn,2:end,indexforDiameter))], 'white'); xticks([1:length(test_i)+3]); xticklabels({'Healthy Data', 'Diabetes Data','No Treatment', 'KN93', 'ML7', 'Y27632', 'CalA', 'CytB'});  ylabel('Fenestration Diameter (nm)')
     B.FaceColor = 'flat';
     B.CData(1,:) = [0 0 1];
     B.CData(2,:) = [0 0 0];
     B.CData(3,:) = [.7 .7 .7];
     for inh = 1:length(test_i)
        B.CData(3+inh,:) = [1 1 1]; 
     end
     hold on;  %er = errorbar([2:7], [mean(s_FC([1:Nn],:,indexforDiameter))], (CI_diameter_lb - CI_diameter_ub)); er.Color = 'r';  er.LineStyle = 'none'; er.LineWidth=2; 
     er = errorbar(1, mean(s_ctrl_d), std(s_ctrl_d)); er.Color = 'r';          er.LineStyle = 'none'; er.LineWidth=2;  
     er = errorbar(2, mean(s_ctrl_db_d), std(s_ctrl_db_d)); er.Color = 'r';          er.LineStyle = 'none'; er.LineWidth=2;
     er = errorbar(3, mean(refForModeld), std(refForModeld)); er.Color = 'r';          er.LineStyle = 'none'; er.LineWidth=2; 
     er = errorbar(4:8, mean(s_FC(1:Nn,2:end,indexforDiameter)),std(s_FC(1:Nn,2:end,indexforDiameter))); er.Color = 'r';  er.LineStyle = 'none'; er.LineWidth=2; %standard deviation instead of CI as the error
     callsigstar(2,p_d_data,'k')
     callsigstar(1,p_c_d,'b')
     ylim([0 120])

     set(gca,'FontSize',8)
    labelstring = {'A', 'B'};
    for v = 1:2
        subplot(2,1,v)
        hold on
        text(-0.1, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize',8)
        set(gca,'FontName','Arial','FontSize',8)
    end

    widthInches = 5.5;
    heightInches = 4.23;
    run('ScriptForExportingImages.m')   

    % [std(s_FC([1:Nn],:,indexforNumber))]
    % [std(s_FC([1:Nn],:,indexforDiameter))]
   
    %% Generate Table E
    % Define chemical agent names
chemical_agents = {'No treatment', 'KN93', 'ML7', 'Y27632', 'Calyculin A', 'Cytochalasin B'};

% Initialize LaTeX table string
latex_table = "\\begin{table}[htbp]\n ";
latex_table = [latex_table, "\\caption{Mean and standard deviation (SD) of predicted fenestration number and diameter after \\textit{in silico} protein knockdown (Fig 6).}\n "];
latex_table = [latex_table, "\\centering\n\\begin{tabular}{lll}\n\\toprule\n "];
latex_table = [latex_table, " & Fenestration Diameter & Fenestration Number \\\\\n "];
latex_table = [latex_table, "Chemical agent & Mean (SD) & Mean (SD) \\\\\n\\midrule\n "];

% Loop through each chemical agent
for i = 1:length(chemical_agents)
    mean_diam = mean(s_FC(:, i, indexforDiameter), 'all');
    std_diam = std(s_FC(:, i, indexforDiameter), 0, 'all');
    mean_num = mean(s_FC(:, i, indexforNumber), 'all');
    std_num = std(s_FC(:, i, indexforNumber), 0, 'all');

    % Format values in scientific notation
    diam_str = sprintf('%.2f ($%.1e$)', mean_diam, std_diam);
    num_str = sprintf('%.3f ($%.1e$)', mean_num, std_num);

    % Append row to LaTeX table
    latex_table = [latex_table, sprintf('%s & %s & %s \\\\\n', chemical_agents{i}, diam_str, num_str)];
end

% Close LaTeX table
latex_table = [latex_table, "\\bottomrule\n\\end{tabular}\n\\label{test_knock_table}\n\\end{table}"];

% Display or write to file
% disp(latex_table);
% Optionally write to file:
fid = fopen('fenestration_table.tex', 'w');
fprintf(fid, '%s', latex_table);
fclose(fid);


     %%
     % figure(6); subplot(2,1,1); b = bar([mean(s_FC([1:Nn],:,indexforNumber))], 'white'); xticks([1:length(test_i)+2]); xticklabels({'healthy', 'KN93', 'ML7', 'Y27632', 'CalA', 'CytB'}); ylabel('Fenestration Number');
     % b.FaceColor = 'flat';
     % b.CData(1,:) = [0 0 1];
     % b.CData(2,:) = [1 1 1];
     % b.CData(3,:) = [1 1 1]; b.CData(4,:) = [1 1 1]; b.CData(5,:) = [1 1 1]; b.CData(6,:) = [1 1 1]; 
     % 
     % set(gca,'FontSize',8)
     % 
     % figure(6); subplot(2,1,2); B = bar([mean(s_FC([1:Nn],:,indexforDiameter))], 'white'); xticks([1:length(test_i)+2]); xticklabels({'healthy',  'KN93', 'ML7', 'Y27632', 'CalA', 'CytB'}); ylabel('Fenestration Diameter (nm)')
     % B.FaceColor = 'flat';
     % B.CData(1,:) = [0 0 1];
     % B.CData(2,:) = [1 1 1];
     % B.CData(3,:) = [1 1 1]; B.CData(4,:) = [1 1 1]; B.CData(5,:) = [1 1 1]; B.CData(6,:) = [1 1 1]; 
     % 
     % set(gca,'FontSize',8)


end
%% 2: LSA-based perturbation
if task == 2

% sensitivity coefficient initialization
% size: length(species), length(params), length(time)
s_FD_Ym = zeros(SP, SP, 1); 
s_FD_W = zeros(SP, RP, 1);

linest = ["--", "-.", ":", "--", "-.", ":",  "--", "-.", ":", "--", "-.", ":", "--", "-.", ":"];

tspan = start_time_h:end_time_h;
percent = 50; opts = [];%odeset(AbsTol=1e-8);


[t, y] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,params,p_params, state, GC_conc', intv); 
y = abs(y);

% Parameter: Ymax
%deltaP = -1; % full knockdown

for m = 1:SP
    params_new = params;
    params_new{3}(m) = 0; %params{3}(m)*(1 + deltaP); % perturb each "Ymax" parameter by a small amount
    opts = [];%odeset(AbsTol=1e-8);
    [time,dy_model] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,params_new,p_params, state, GC_conc', intv);
    

    dy_model = real(dy_model);

    dy_modelR = real(dy_model);
    
    
    for l = 1:size(dy_modelR,2)
        s_FD_Ym(l,m,:) = (dy_modelR(end,l) - y(end,l))/(params_new{3}(m) - params{3}(m))* params{3}(m)/y(end,l)*100; %/(percent*1e-2); % at 20 weeks
    end
%     figure(2);
%     subplot(5,7,m);
%     plot(time/(24*7), dy_modelR(:,indexforNumber)); hold on; legend(params{4}(m))
end


%%
% Parameter: W
for m = 1:RP
    params_new = params;
    params_new{1}(1,m) = 0; %params{1}(1,m)*(1 + deltaP); % perturb each "W" parameter by a small amount
        
    [time,dy_model] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,params_new,p_params, state, GC_conc', intv);
     
    dy_model = real(dy_model);
    
    
    dy_modelR = real(dy_model);
    
    for l = 1:size(dy_modelR,2)
        s_FD_W(l,m,:) = ((dy_modelR(end,l)) - y(end,l))/(params_new{1}(1,m) - params{1}(1,m))* params{1}(1,m)/y(end,l)*100; %/(percent*1e-2); % difference in values at 20 weeks
    end
end

figure(55); 
hcolormap = colMapGen([1 0 0],[0 0 1],100,1);
% do not consider the GLU input i = 1 or the responses Number and Diameter i = 36 and 37 because species parameters are undefined
heatmapdata = real(s_FD_Ym(indexforNumber,2:35,1));
[x_sorted, sortIdx] = sort(heatmapdata);
subplot(2,1,1); h1=heatmap(x_sorted, 'Colormap', hcolormap); 
ax = gca; ax.Interpreter = 'tex';
% sortIdx is offset by 1 because heatmap started at 2 to ignore GLU
ax.XDisplayLabels = params{4}(sortIdx+1); ax.XLabel = 'Species i';
ax.YDisplayLabels = 'y_{max_i}'; ax.YLabel = {'Normalized % change'; 'in Number'; 'relative to'};
%params{4}(indexforNumber); 
set(gca,'FontName','Arial','FontSize',6)
h1.ColorLimits = [-round(max(abs(x_sorted)),-1), round(max(abs(x_sorted)),-1)];
h1.Title = 'A';
filtered_sorted_indices = find(abs(x_sorted) > 5);
original_indices = sortIdx(filtered_sorted_indices);
numSpeciesSortIdx = original_indices+1;

hcolormap = colMapGen([1 0 0],[0 0 1],50,1);
% do not consider the GLU input rxn j = 1 because the reaction parameters are undefined
heatmapdata = real(s_FD_W(indexforNumber,2:47,1));
[x_sorted, sortIdx] = sort(heatmapdata);
subplot(2,1,2); h2=heatmap(x_sorted, 'Colormap', hcolormap); 
ax = gca; ax.Interpreter = 'tex';
% sortIdx is offset by 1 because heatmap started at 2 to ignore GLU
ax.XDisplayLabels = params{5}(sortIdx+1); ax.XLabel = 'Reaction rules j';
ax.YDisplayLabels = 'W_j'; ax.YLabel = {'Normalized % change'; 'in Number'; 'relative to'};
%params{4}(indexforNumber);  
set(gca,'FontName','Arial','FontSize',6)
h2.Title = 'B';
h2.ColorLimits = [-round(max(abs(x_sorted)),-1), round(max(abs(x_sorted)),-1)];
filtered_sorted_indices = find(abs(x_sorted) > 5);
original_indices = sortIdx(filtered_sorted_indices);
numRxnSortIdx = original_indices+1;

hcolormap = colMapGen([1 0 0],[0 0 1],80,0.5);
% do not consider the GLU input i = 1 or the responses Number and Diameter i = 36 and 37 because species parameters are undefined
heatmapdata = real(s_FD_Ym(indexforDiameter,2:35,1)); 
[x_sorted, sortIdx] = sort(heatmapdata);
figure(56); subplot(2,1,1); h1=heatmap(x_sorted, 'Colormap', hcolormap); 
ax = gca; ax.Interpreter = 'tex';
% sortIdx is offset by 1 because heatmap started at 2 to ignore GLU
ax.XDisplayLabels = params{4}(sortIdx+1); ax.XLabel = 'Species i';
ax.YDisplayLabels = 'y_{max_i}'; ax.YLabel = {'Normalized % change'; 'in Diameter'; 'relative to'};
%params{4}(indexforDiameter)  ; 
set(gca,'FontName','Arial','FontSize',6)
h1.Title = 'A';
h1.ColorLimits = [0, round(max(abs(x_sorted)),-1)];
filtered_sorted_indices = find(abs(x_sorted) > 5);
original_indices = sortIdx(filtered_sorted_indices);
diamSpeciesSortIdx = original_indices+1;

hcolormap = colMapGen([1 0 0],[0 0 1],80,0.5);
% do not consider the GLU input rxn j = 1 because the reaction parameters are undefined
heatmapdata = real(s_FD_W(indexforDiameter,2:47,1));
[x_sorted, sortIdx] = sort(heatmapdata);
subplot(2,1,2); h2=heatmap(x_sorted, 'Colormap', hcolormap); 
ax = gca; ax.Interpreter = 'tex';
% sortIdx is offset by 1 because heatmap started at 2 to ignore GLU
ax.XDisplayLabels = params{5}(sortIdx+1);  ax.XLabel = 'Reaction rules j';
ax.YDisplayLabels = 'W_j'; ax.YLabel = {'Normalized % change'; 'in Diameter'; 'relative to'};
%params{4}(indexforDiameter); 
set(gca,'FontName','Arial','FontSize',6)
h2.Title = 'B';
h2.ColorLimits = [0, round(max(abs(x_sorted)),-1)];
filtered_sorted_indices = find(abs(x_sorted) > 5);
original_indices = sortIdx(filtered_sorted_indices);
diamRxnSortIdx = original_indices+1;

figure(55)
figname = 'FigE';
widthInches = 5;
heightInches = 4.23;
run('ScriptForExportingImages.m')     

figure(56)
figname = 'FigF';
widthInches = 5;
heightInches = 4.23;
run('ScriptForExportingImages.m')

listOfSensSortIdx = {numSpeciesSortIdx,numRxnSortIdx,diamSpeciesSortIdx,diamRxnSortIdx};
save('data/listOfSensSortIdx.mat','listOfSensSortIdx')

end

%%  3: time-dependent intervention
if task == 3

    % sensitivity coefficient initialization
    % size: length(species), length(params), length(time)
    s_FD_Ym = zeros(SP, SP, 1); 
    s_FD_W = zeros(SP, RP, 1);
    
    linest = ["--", "-.", ":", "--", "-.", ":",  "--", "-.", ":", "--", "-.", ":", "--", "-.", ":","--", "-.", ":"];
    
    tspan = start_time_h:1:Tstop;
    percent = 50; 
    opts = odeset(AbsTol=1e-8,RelTol=1e-6);
    % In this task we use ode23s and stricter than default ode tolerances to provide more stability when the
    % perturbation cause discontinuities in parameters at differen times
    longTfinal = 30; %weeks
    longTfinal_h = longTfinal*24*7;
    [t, y] = ode23s(@coupledODE_IVV_step,tspan,y0,opts,params,p_params, state, GC_conc', intv); 
    y = abs(y);

    [tfull, yfull] = ode23s(@coupledODE_IVV_step,[start_time_h:longTfinal_h],y0,opts,params,p_params, state, GC_conc', intv); 
    yfull = abs(yfull);
    
    if Tstop/24/7 == 8
        NUM = 7; % FIG 7 publication
        figname = 'Fig7';
    elseif Tstop/24/7 == 10
        NUM = 57; % FIG G publication
        figname = 'FigG';
    elseif Tstop/24/7 == 20
        NUM = 58; % FIG H publication
        figname = 'FigH';
    end
load('data/listOfSensSortIdx.mat','listOfSensSortIdx');
% {numSpeciesSortIdx,numRxnSortIdx,diamSpeciesSortIdx,diamRxnSortIdx};
    figure(NUM); hold on; ax = gca; ax.FontSize = 8;
    subplot(2,2,1); box; hold on; h=plot(tfull/(24*7), yfull(:,indexforNumber), 'LineWidth', 1.2, 'color', 'k','DisplayName','No inhibition'); 
    legendHandles = h;
    legendNames = {'No inhibition'};
   
    inhib_knock_n = listOfSensSortIdx{1}; %all: 2:35; %KP list: [1,2,3,5,7,9,11,12,18,19,25,27,31,33]; % 6 ignored
    inhib_prod_n = listOfSensSortIdx{2}; %all: 2:47; %KP list: [46,21,24,16,2,3,18,20,29,41,43];  % 17 ignored
    % % inhib_knock_n2 = [35,36,31]; %promotes fenestrations
    
    tnew = Tstop:1:longTfinal_h;
    y0(:) = y(end,:);
    z_params = params;
    for inh = 1:length(inhib_knock_n)
        z_params{3}(inhib_knock_n(inh)) = 0.5;
        [tn, yn] = ode23s(@coupledODE_IVV_step,tnew,y0,opts,z_params,p_params, state, GC_conc', intv); 
        yn = abs(yn);
        
       
        name = params{4}(inhib_knock_n(inh));
        h=plot(tn/(24*7), yn(:,indexforNumber), 'LineWidth', 1.2,'LineStyle',linest(inh), 'DisplayName',name{1}); xlabel('Time (week)'); ylabel('Fenestration Number');
        legendHandles(end+1) = h;
        legendNames{end+1} = name{1};

        z_params = params;
    end

    legend(legendHandles, legendNames, 'Location', 'Eastoutside','fontsize',4)

    subplot(2,2,2); box; hold on; h=plot(tfull/(24*7), yfull(:,indexforNumber), 'LineWidth', 1.2, 'color', 'k','DisplayName','No inhibition'); 
    legendHandles = h;
    legendNames = {'No inhibition'};
    z_params = params;
    for inh = 1:length(inhib_prod_n)
        %z_params{3}(inhib_knock_n2(inh)) = 2;
        z_params{1}(1,inhib_prod_n(inh)) = 0.5*z_params{1}(1,inhib_prod_n(inh));
        [tn, yn] = ode23s(@coupledODE_IVV_step,tnew,y0,opts,z_params,p_params, state, GC_conc', intv); 
        yn = abs(yn);
        
        
        name = params{5}(inhib_prod_n(inh));
        h=plot(tn/(24*7), yn(:,indexforNumber), 'LineWidth', 1.2,'LineStyle',linest(inh), 'DisplayName',name{1}); xlabel('Time (week)'); ylabel('Fenestration Number');
        legendHandles(end+1) = h;
        legendNames{end+1} = name{1};        
        z_params = params;
    end
    legend(legendHandles, legendNames, 'Location', 'Eastoutside','fontsize',4)

    
    inhib_knock_d = listOfSensSortIdx{3}; % all: 2:35;%KP list: [1,3,5,6,7,9,11,12,18,19,25,27,31,33]; 
    inhib_prod = listOfSensSortIdx{4}; % all: 2:47;%KP list: [2,3,19,20,29,41,43,15,24,21];
   
    
    subplot(2,2,3); box; hold on; h=plot(tfull/(24*7), yfull(:,indexforDiameter), 'LineWidth', 1.2, 'color', 'k','DisplayName','No inhibition'); 
    legendHandles = h;
    legendNames = {'No inhibition'};
    z_params = params;
    for inh = 1:length(inhib_knock_d)
        z_params{3}(inhib_knock_d(inh)) = 0.5;
        [tn, yn] = ode23s(@coupledODE_IVV_step,tnew,y0,opts,z_params,p_params, state, GC_conc', intv); 
        yn = abs(yn);
              
        name = params{4}(inhib_knock_d(inh));
        h=plot(tn/(24*7), yn(:,indexforDiameter), 'LineWidth', 1.2,'LineStyle',linest(inh), 'DisplayName',name{1}); xlabel('Time (week)'); ylabel('Fenestration Diameter (nm)');
        legendHandles(end+1) = h;
        legendNames{end+1} = name{1};
        z_params = params;
    end
    
    legend(legendHandles, legendNames, 'Location', 'Eastoutside','fontsize',4)

    subplot(2,2,4); box; hold on; h=plot(tfull/(24*7), yfull(:,indexforDiameter), 'LineWidth', 1.2, 'color', 'k','DisplayName','No inhibition');  
    legendHandles = h;
    legendNames = {'No inhibition'};


    z_params = params;
    for inh = 1:length(inhib_prod)
        z_params{1}(1,inhib_prod(inh)) = 0.5*z_params{1}(1,inhib_prod(inh));
        [tn, yn] = ode23s(@coupledODE_IVV_step,tnew,y0,opts,z_params,p_params, state, GC_conc', intv); 
        yn = abs(yn);
        
        
        name = params{5}(inhib_prod(inh));
        h=plot(tn/(24*7), yn(:,indexforDiameter), 'LineWidth', 1.2,'LineStyle',linest(inh), 'DisplayName',name{1}); xlabel('Time (week)'); ylabel('Fenestration Diameter (nm)');
        legendHandles(end+1) = h;
        legendNames{end+1} = name{1};
        z_params = params;
    end

    legend(legendHandles, legendNames, 'Location', 'Eastoutside','fontsize',4)


    labelstring = {'A', 'B','C','D'};
    for v = 1:4
        subplot(2,2,v)
        hold on
        text(-0.1, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize',8)
        set(gca,'FontName','Arial','FontSize',8)
    end

    widthInches = 9;
    heightInches = 5.5;
    run('ScriptForExportingImages.m')   

end
%%  4: Glucose-intervention
if task == 4
    opts = [];%odeset(RelTol=1e-6);
    % no intervention step
    intv = 'none';
    t1 = start_time_h:end_time_h;
    %[T, Y] = ode15s(@coupledODE_IVV_step,t1,y0,opts,params,p_params, state, GC_conc', intv);
    %Y = real(Y);

    Nn = 100;
    glu_sampled = zeros(11,Nn);
    rng("twister") % Default random number generator algorithm with seed = 0 to ensure that we generate the same sequence of draws
    glucose_data_Lee_sd = abs(glucose_lee - LB_lee);
    glucose_data_Finch_sd = abs(glu_finch - glu_UB); 
    subsetIdx = find(time_lee>=6*24*7);
    sigma_data = [glucose_data_Lee_sd(subsetIdx)'  glucose_data_Finch_sd];

%     for Nstep = 1:Nn
%         for i = 1:length(GC_time)
% %             glu_sampled(i,Nstep) = unifrnd(GC_LB(:,i), GC_UB(:,i)); %            
%             glu_sampled(i,Nstep) = normrnd(GC_conc(:,i), sigma_data(i)); %
%         end
    for i = 1:length(GC_time)
%                glu_sampled(i,Nstep) = unifrnd(GC_LB(:,i), GC_UB(:,i)); %
        glu_sampled(i,:) = normrnd(GC_conc(i), sigma_data(i),[1,Nn]); %
    end
    for Nstep = 1:Nn

        [T, Y] = ode15s(@coupledODE_IVV_step,t1,y0,opts,params,p_params, state, glu_sampled(:,Nstep), intv);
        YstepP(Nstep,:,:) = real(Y);
        %fprintf('run %i finished\n', Nstep)

    end

    Ymean(1,:,:) = mean(YstepP(1:Nn,:,:));
    Ymean = squeeze(Ymean(1,:,:));   
       

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

    % intervention at 4 hours 
    intv = '4h';
    [T4, Y4] = ode15s(@coupledODE_IVV_step,t1,y0,opts,params,p_params, state, GC_conc', intv);
    Y4 = real(Y4);

    
       % Gp = step_function(glu_sampled);
        for Tt = start_time_h:1:tspan(end)
            if (Tt >= time_lee(1) && Tt <= time_lee(3))
                GLU_p4(Tt,1) = 0.051*(Tt)  - 9.38;
            else
                GLU_p4(Tt,1) = 0.051*(start_time_h)  - 9.38;
            end
        end
    
    y0 = Y4(end,:);
    intv = '10h';
    % intervention at 10 hours 
    [T10, Y10] = ode23s(@coupledODE_IVV_step,t1,y0,opts,params,p_params, state, GC_conc', intv);
    Y10 = real(Y10);

    opts = odeset(AbsTol=1e-8);
    Nn ;
    glu_sampled = zeros(11,Nn);
    rng("twister") % Default random number generator algorithm with seed = 0 to ensure that we generate the same sequence of draws
    glucose_data_Lee_sd = abs(glucose_lee - LB_lee);
    glucose_data_Finch_sd = abs(glu_finch - glu_UB); 
    subsetIdx = find(time_lee>=6*24*7);
    sigma_data = [glucose_data_Lee_sd(subsetIdx)'  glucose_data_Finch_sd];

%     for Nstep = 1:Nn
%         for i = 1:length(GC_time)
% %             glu_sampled(i,Nstep) = unifrnd(GC_LB(:,i), GC_UB(:,i)); %            
%             glu_sampled(i,Nstep) = normrnd(GC_conc(:,i), sigma_data(i)); %
%         end
    for i = 1:length(GC_time)
%                glu_sampled(i,Nstep) = unifrnd(GC_LB(:,i), GC_UB(:,i)); %
        glu_sampled(i,:) = normrnd(GC_conc(i), sigma_data(i),[1,Nn]); %
    end
    for Nstep = 1:Nn

        [T10, Y10] = ode23s(@coupledODE_IVV_step,t1,y0,opts,params,p_params, state, glu_sampled(:,Nstep), intv);
        YstepP10(Nstep,:,:) = real(Y10);
        %fprintf('run %i finished\n', Nstep)

    end

    Ymean10(1,:,:) = mean(YstepP10(1:Nn,:,:));
    Ymean10 = squeeze(Ymean10(1,:,:));   
    
    for Nstep=1:Nn
        %Gp = step_function(glu_sampled(:,Nstep));
        for Tt = start_time_h:1:tspan(end)
            if (Tt >= time_lee(1) && Tt <= time_lee(5))
                GLU_p10(Tt,Nstep) =  0.051*(Tt)  - 9.38;
            elseif (Tt>time_lee(5) && Tt<=time_lee(9))
           
                GLU_p10(Tt,Nstep) = step_function(Tt, glu_sampled(:,Nstep));
            else
                GLU_p10(Tt,Nstep) = 0.051*(start_time_h)  - 9.38;
            end
        end
    end

    
 %%   

    figure(4); 
    figname = 'Fig4';

    Gp0 = GLU_p(start_time_h,1);
    leftymin = 4;
    leftymax = 55;
    % need to convert the y axis to normalized units to get the scaling
    % perfect
    rightymin = (leftymin - Gp0) / (max(glu_UB) - Gp0);
    rightymax = (leftymax - Gp0) / (max(glu_UB) - Gp0);

    box;  
    
    % subplot(1,2,1)
    % yyaxis left
    % hold on
    % 
    % set(gca,'FontName','Arial','FontSize',8)
    % 
    % % colors = jet(25);
    % % colororder(colors);
    % colororder('default');
    % 
    % for i = 1:25;%size(GLU_p,2)
    %     plot(T/(24*7), GLU_p(start_time_h:end_time_h,i), 'LineWidth', 1.2,'LineStyle','-','DisplayName',num2str(i));
    % end
    % % % for aesthetics make the mean black
    % % MEANGLU_p=mean(GLU_p,2);
    % % h=plot(T/(24*7),  MEANGLU_p, 'LineWidth', 5,'LineStyle','-','Color','k');
    % % %legend('show', 'NumColumns', 3, 'Location', 'southoutside','FontSize',4);
    % % h.DisplayName=sprintf('Mean of samples = \nModel: no intervention case')
    % % legend('show', 'NumColumns', 8,'Location', 'southoutside','FontSize',4);
    % 
    % 
    % ax = gca;
    % ax.YAxis(1).Color = 'k'; 
    % hold on;
    % xlabel('Time (weeks)'); xlim([0,21]);
    % ylabel('Glucose (mmol/l)'); 
    % 
    % ylim([leftymin,leftymax]);
    % 
    % yyaxis right
    % set(gca,'FontName','Arial','FontSize',8)
    % ylabel('Normalized Glucose Units');
    % ax.YAxis(2).Color = 'k';    
    % ylim([rightymin,rightymax]);

    % subplot(1,2,2)
    yyaxis left
    hold on
    plot(T/(24*7), mean(GLU_p(start_time_h:end_time_h,:)'), 'LineWidth', 1.2, 'Color', 'k'); 
    plot(T4/(24*7), GLU_p4(start_time_h:end_time_h,1), 'LineWidth', 1.2, 'Color', 'k', 'LineStyle', ':'); 
    plot(T10/(24*7), mean(GLU_p10(start_time_h:end_time_h,:)'), 'LineWidth', 1.2, 'Color', [0.7 0.7 0.7], 'LineStyle', '--'); 
    glucose_data_Lee_sd = abs(glucose_lee - LB_lee);
    glucose_data_Finch_sd = abs(glu_finch - glu_UB); 
    hold on; errorbar(time_lee/(7*24), glucose_lee, glucose_data_Lee_sd, 'o', 'Color',[0  0  1],'MarkerFaceColor',[0  0  1]);
    hold on; errorbar(time_finch/(7*24), glu_finch, glucose_data_Finch_sd, '^', 'Color',[1 0 0],'MarkerFaceColor',[1 0 0]);

    set(gca,'FontName','Arial','FontSize',8)
    ax = gca;
    ax.YAxis(1).Color = 'k'; 
    hold on;
    xlabel('Time (weeks)'); xlim([0,21]);
    ylabel('Glucose (mmol/l)'); 
    
    ylim([leftymin,leftymax]);

    yyaxis right
    set(gca,'FontName','Arial','FontSize',8)
    ylabel('Normalized Glucose Units');
    ax.YAxis(2).Color = 'k';    
    ylim([rightymin,rightymax]);
    box on

    % labelstring = {'A', 'B'};
    % for v = 1:2
    %     subplot(1,2,v)
    %     hold on
    %     text(-0.15, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize',8)
    %     set(gca,'FontName','Arial','FontSize',8)
    % end

    legend('No intervention', 'Intervention at 4 weeks', 'Intervention at 10 weeks', 'Lee et al. (2018)','Finch et al. (2022)','location','northeast');
    widthInches = 5.5;
    heightInches = 4.23;
    run('ScriptForExportingImages.m')    


%     subplot(1,2,2); 
%     plot(T/(24*7), Y(:,1), 'LineWidth', 3, 'Color', 'k');
%     hold on;
%     plot(T4/(24*7), Y4(:,1), 'LineWidth', 3,'Color', 'k', 'LineStyle', ':'); 
%     hold on;
%     plot(T10/(24*7), Y10(:,1), 'LineWidth', 3, 'Color', [0.7 0.7 0.7], 'LineStyle', '--'); 
%     xlabel('Time (weeks)');
%     ylabel('Normalized Glucose');% xlim([0,22]); 
%     ax = gca; ax.FontSize = 20;
%     hold on
%     x0=10;
%     y0=10;
%     width=800;
%     height=800;
%     set(gcf,'position',[x0,y0,width,height])
%     legend('Model')
%     box on
    

%%
figure(5);
figname = 'Fig5';
    subplot(3,3,[1 2]);
    bar((Ymean(end,2:35) - Ymean(1,2:35)), 'k'); title('No intervention'); 
    xticks(1:34); xticklabels(params{4}(2:35)); ylabel('Change relative to baseline')
    
    subplot(3,3,3)
    bar((Ymean(end,indexforNumber:indexforDiameter) - Ymean(1,indexforNumber:indexforDiameter)), 'k'); title('No intervention')
    xticks(1:2); xticklabels(params{4}(indexforNumber:indexforDiameter)); ylabel('Change relative to baseline')
 
    subplot(3,3,[4,5]);
    bar((Y4(end,2:35) - Y4(1,2:35)), 'k'); title('Intervention at 4 weeks'); 
    xticks(1:34); xticklabels(params{4}(2:35)); ylabel('Change relative to baseline')
    
    subplot(3,3,6)
    bar((Y4(end,indexforNumber:indexforDiameter) - Y4(1,indexforNumber:indexforDiameter)), 'k'); title('Intervention at 4 weeks')
    xticks(1:2); xticklabels(params{4}(indexforNumber:indexforDiameter)); ylabel('Change relative to baseline')
    
    subplot(3,3,[7,8]);
    bar((Ymean10(end,2:35) - Ymean10(1,2:35)), 'k'); title('Intervention at 10 weeks'); 
    xticks(1:34); xticklabels(params{4}(2:35)); ylabel('Change relative to baseline')
    
    subplot(3,3,9)
    bar((Ymean10(end,indexforNumber:indexforDiameter) - Ymean10(1,indexforNumber:indexforDiameter)), 'k'); title('Intervention at 10 weeks')
    xticks(1:2); xticklabels(params{4}(indexforNumber:indexforDiameter)); ylabel('Change relative to baseline')

    labelstring = {'A', 'B','C', 'D','E', 'F'};
    label_iterator = 0;
    for v = 1:9
        if mod(v,3) == 1
            subplot(3,3,[v v+1]);
            hold on
        
            label_iterator = label_iterator+1;
            text(-0.1, 1.1, labelstring(label_iterator)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize',8)
            set(gca,'FontName','Arial','FontSize',4)
        elseif mod(v,3) == 0
            subplot(3,3,v)
            hold on
            label_iterator = label_iterator+1;
            text(-0.1, 1.1, labelstring(label_iterator)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize',8)
            set(gca,'FontName','Arial','FontSize',4)      
        end
    end

    widthInches = 5.5;
    heightInches = 5.5;
    run('ScriptForExportingImages.m')   
end



end
   function output = callsigstar(ctrlgroupBarNumber,pvalues,color)
       % call sigstart for simulation groups that are significant based on
       % input p-value but don't create overbars for non-significant values

       % don't do self comparisons or comparisons to the left that have already been done
       pvalues = pvalues(ctrlgroupBarNumber+1:end);
       pvalues < 0.05
       pvalues < 0.01
       pvalues < 0.001
       pvalues < 0.0001
       significant_indices = find(pvalues < 0.05);

       Groups = arrayfun(@(x) [ctrlgroupBarNumber, x+ctrlgroupBarNumber], significant_indices, 'UniformOutput', false)

       H = sigstar(Groups,pvalues(significant_indices),1);
       set(H,'color',color)
   end