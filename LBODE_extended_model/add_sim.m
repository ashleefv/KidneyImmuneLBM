% Additional Simulation using Healthy Mice Glucose Data from Finch et al. 
% Changes in mean fenestration number and diameter are plotted as barplot
% before and after treatment with agent

function [s_FC] = add_sim(params, y0, tspan, p_params, state, task)
global number_ctrl time_ctrl density diameter GC_conc GC_time GC_LB GC_UB time_lee glu_UB glu_LB time_finch glu_finch LB_lee UB_lee glucose_lee ctrl_glu ctrl_LB ctrl_UB
G = load('data/GLU_data.mat');                          % Glucose concentration from Finch et al. 2022 (12-20 weeks) and Lee et al. 2018 (2-11 weeks). Fig 1b male ob-/ob- from both sources


% Time (long-term mice sim.)
    start_time = 2; %weeks
    start_time_h = start_time*7*24;
    end_time = 20; %weeks
    end_time_h = end_time*7*24;
    tspan = start_time_h:1:end_time_h; % hours
indexforNumber = 36;
indexforDiameter = 37;

intv = 'none';

z_params = params;
%% 1: test_knockout
if task == 1
    opts=[];
    Nn = 100;
    state = state;
    test_i = [29,32,31,34,2]; % inhibitors: KN93, ML7, Y27632, CalA, CytB

    rng("twister") % Default random number generator algorithm with seed = 0 to ensure that we generate the same sequence of draws


    if state == "norm_mice"
        
       % glu_sampled = zeros(11,1);
       % glu_sampled([7:11],1) = ctrl_glu;

        glu_sampled = zeros(11,Nn);
        glu_ctrl_ref=readmatrix('data/LEE_FINCH_CTRL_GLU.csv');
        mean_ctrl = glu_ctrl_ref([5:end],2);
        sd_ctrl = glu_ctrl_ref([5:end],4);

        for i = 1:length(GC_time)
            glu_sampled(i,:) = normrnd(mean_ctrl(i), sd_ctrl(i), [1,Nn]); %
        end

        for Nstep = 1:Nn
        
            [t, y] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,params,p_params, state, glu_sampled(:,Nstep), intv);
            YstepP(Nstep,:) = real(y(end,:));
            %disp(size(YstepP))
        end
    end
    s_FC(:,1,:) = YstepP(:,:); % no treatment case

 %% read data sets for healthy and diabetes at 20 weeks with the samples from Finch
    % FINCH figure 1 E & F
    s_ref=readmatrix('data/FINCH_FENESTRATION_20wk.csv');
    s_ctrl_n = s_ref(:,2); % number data for healthy case
    s_ctrl_d = s_ref(:,3); % diameter data for healthy case
    % s_ctrl_db_n = s_ref(1:8,4); % number data for diabetic case
    % s_ctrl_db_d = s_ref(1:8,5); % diameter data for diabetic case
    % s_ctrl_n([1:Nn],1,1) = 5.8; % reference is set for initial data value for healthy case
    % s_ctrl_d([1:Nn],1,1) = 50;  % reference is set for initial data value for healthy case
    test_i = [29,32,31,34,2]; % inhibitors: KN93, ML7, Y27632, CalA, CytB
    z_params = params;
    % we plot 3 additional columns before test_i: healthy data, diabetes
    % data, no treatment

    % % These are the posterior distributions of the output fenestration number and diameter based on MC sampling of the parameter posteriors    
    % MC = load('data/MC_25_fen_multirun.mat');
    % Y_param_var = MC.Y_param_var;
    % MC_samplesn = Y_param_var(:,end,indexforNumber+1);
    % MC_samplesd = Y_param_var(:,end,indexforDiameter+1);

    % run t-tests for comparison of number and diameter to each of these
    % three columns
    testColumns = 1;
    p_c_n = zeros(1,length(test_i)+testColumns); 
    p_c_d = zeros(1,length(test_i)+testColumns);
    p_n_data = zeros(1,length(test_i)+testColumns);
    p_d_data = zeros(1,length(test_i)+testColumns);
    p_n = zeros(1,length(test_i)+testColumns);
    p_d = zeros(1,length(test_i)+testColumns);
    offset = 1;

    %%

    % We have two options for a distribution of samples of simulated untreated
% fenestration number and diameter. One of the MC_samples* were * is n or d
% for number and diameter. The second option selected here is the no
% treatment case sampled over glucose inputs. The latter is closer to the
% normal distribution that is expected in a t-test, so is used here.
refForModeln =  s_FC(:,1,indexforNumber);
refForModeld =  s_FC(:,1,indexforDiameter);

    % % compare diabetes data to control data and store in the offset column
    % %p-values for ctrl number relative to the healty condition no treatment simulation
    % [h, p_c_n(1,offset)] = ttest2(s_ctrl_db_n, s_ctrl_n, 'Alpha', 0.05,'Vartype','unequal');
    % % %p-values for ctrl diameter relative to the healty condition no treatment
    % [h, p_c_d(1,offset)] = ttest2(s_ctrl_db_d, s_ctrl_d, 'Alpha', 0.05,'Vartype','unequal');
    

    inh = 0; % no treatment case
        %p-values for ctrl number relative to the healty condition no treatment data
        [h, p_c_n(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforNumber), s_ctrl_n, 'Alpha', 0.05,'Vartype','unequal');
        %p-values for ctrl diameter relative to the healty condition no treatment data
        [h, p_c_d(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforDiameter), s_ctrl_d, 'Alpha', 0.05,'Vartype','unequal');

        
        %  %p-values for number relative to the diseased condition no treatment DATA
        % [h, p_n_data(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforNumber), s_ctrl_db_n, 'Alpha', 0.05,'Vartype','unequal');
        % %p-values for diameter relative to the diseased condition no treatment DATA
        % [h, p_d_data(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforDiameter), s_ctrl_db_d, 'Alpha', 0.05,'Vartype','unequal');

        %[h, p_n(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforNumber), refForModeln, 'Alpha', 0.05,'Vartype','unequal');
        %p-values for diameter relative to the diseased condition no treatment simulation
        %[h, p_d(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforDiameter), refForModeld, 'Alpha', 0.05,'Vartype','unequal');
        
       

    for inh = 1:length(test_i)
        % treatment knocks out a pathway via its parameter
        z_params{3}(test_i(inh)) = 0;
        
       % glu_sampled = zeros(11,1);
       % glu_sampled([7:11],1) = ctrl_glu;


        for Nstep = 1:Nn


            
            [tn, yn] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,z_params,p_params, state,  glu_sampled(:,Nstep), intv);
            YstepP_new(Nstep,:) = real(yn(end,:));
            %fprintf('run %i finished\n', Nstep)

        end
        fprintf('finished treatment %i \n', inh)
       
       
        s_FC(:,inh+1,:) = YstepP_new(:,:);


        

        %p-values for  number relative to the healty condition no treatment data
        [h, p_c_n(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforNumber), s_ctrl_n, 'Alpha', 0.05,'Vartype','unequal');
        %p-values for diameter relative to the healty condition no treatment data
        [h, p_c_d(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforDiameter), s_ctrl_d, 'Alpha', 0.05,'Vartype','unequal');
        
        %p-values for number relative to the diseased condition no treatment DATA
        % [h, p_n_data(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforNumber), s_ctrl_db_n, 'Alpha', 0.05,'Vartype','unequal');
        % %p-values for diameter relative to the diseased condition no treatment DATA
        % [h, p_d_data(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforDiameter), s_ctrl_db_d, 'Alpha', 0.05,'Vartype','unequal');
        %[h, p_n(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforNumber), refForModeln, 'Alpha', 0.05,'Vartype','unequal');
        %p-values for diameter relative to the diseased condition no treatment simulation
        %[h, p_d(1,offset+inh+1)] = ttest2(s_FC(:,inh+1,indexforDiameter), refForModeld, 'Alpha', 0.05,'Vartype','unequal');
        

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

%%
    s_FC_d = squeeze(s_FC(:,:,indexforDiameter));
    s_FC_n = squeeze(s_FC(:,:,indexforNumber));

    CI_diameter_ub = mean(s_FC_d,1) + std(s_FC_d,1)*1.96; % upper-bound for diameter predicted
    CI_diameter_lb = mean(s_FC_d,1) - std(s_FC_d,1)*1.96;

    CI_number_ub = mean(s_FC_n,1) + std(s_FC_n,1)*1.96; % upper-bound for number predicted
    CI_number_lb = mean(s_FC_n,1) - std(s_FC_n,1)*1.96; 
%%
%     disp(mean(s_FC_d,1)); disp(std(s_FC_d,1));
%     disp(mean(s_FC_n,1)); disp(std(s_FC_n,1));
% 
% %%
%      s_FC([1:Nn],:,indexforNumber)
%      size(s_FC)
     figure(6)
     figname = 'Fig6norm';
     subplot(2,1,1); b = bar([mean(s_ctrl_n), mean(refForModeln), mean(s_FC(1:Nn,2:end,indexforNumber))], 'white'); xticks(1:length(test_i)+2); 
    xticklabels({'Healthy Data', 'Healthy Model', 'KN93', 'ML7', 'Y27632', 'CalA', 'CytB'}); 
     ylabel('Fenestration Number');
     b.FaceColor = 'flat';
     b.CData(1,:) = [0 0 1];
     b.CData(2,:) = [.7 .7 .7];
     for inh = 1:length(test_i)
        b.CData(2+inh,:) = [1 1 1]; 
     end
   hold on;  
   %er = errorbar([2:7], [mean(s_FC([1:Nn],:,indexforNumber))], (CI_number_lb-CI_number_ub)); er.Color = 'r';          er.LineStyle = 'none'; er.LineWidth=2;  
  er = errorbar(1, mean(s_ctrl_n), std(s_ctrl_n)); er.Color = 'r';          er.LineStyle = 'none'; er.LineWidth=2; 
  % er = errorbar(2, mean(refForModeln), std(refForModeln)); er.Color = 'r';          er.LineStyle = 'none'; er.LineWidth=2;    
  % er = errorbar(3:7, mean(s_FC(1:Nn,2:end,indexforNumber)), std(s_FC(1:Nn,2:end,indexforNumber))); er.Color = 'r';          er.LineStyle = 'none'; er.LineWidth=2;  
   ylim([0 10.5])
   
    % callsigstar(1,p_c_n,'b')

     set(gca,'FontSize',8)

     figure(6); subplot(2,1,2); B = bar([mean(s_ctrl_d), mean(refForModeld), mean(s_FC(1:Nn,2:end,indexforDiameter))], 'white'); xticks([1:length(test_i)+2]); 
     xticklabels({'Healthy Data', 'Healthy Model', 'KN93', 'ML7', 'Y27632', 'CalA', 'CytB'});  ylabel('Fenestration Diameter (nm)')
          B.FaceColor = 'flat';
     B.CData(1,:) = [0 0 1];
     B.CData(2,:) = [.7 .7 .7];
     for inh = 1:length(test_i)
        B.CData(2+inh,:) = [1 1 1]; 
     end
     hold on;  %er = errorbar([2:7], [mean(s_FC([1:Nn],:,indexforDiameter))], (CI_diameter_lb - CI_diameter_ub)); er.Color = 'r';  er.LineStyle = 'none'; er.LineWidth=2; 
    er = errorbar(1, mean(s_ctrl_d), std(s_ctrl_d)); er.Color = 'r';          er.LineStyle = 'none'; er.LineWidth=2;  
    % er = errorbar(2, mean(refForModeld), std(refForModeld)); er.Color = 'r';          er.LineStyle = 'none'; er.LineWidth=2; 
    % er = errorbar(3:7, mean(s_FC(1:Nn,2:end,indexforDiameter)),std(s_FC(1:Nn,2:end,indexforDiameter))); er.Color = 'r';  er.LineStyle = 'none'; er.LineWidth=2; %standard deviation instead of CI as the error
    % 
    % callsigstar(1,p_c_d,'b')
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
 %% HEALTHY MICE SIM
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