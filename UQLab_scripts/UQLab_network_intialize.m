%% SENSITIVITY: COMPARISON OF MC- VS. PCE- AND LRA-BASED SOBOL' INDICES
%
% In this example, Sobol' sensitivity indices for the 
% 
% are calculated with three different methods:
% Monte Carlo (MC) simulation, polynomial chaos expansion (PCE), 
% and canonical low-rank approximation (LRA).

%% 1 - INITIALIZE UQLAB
%
% Clear variables from the workspace,
% set random number generator for reproducible results,
% and initialize the UQLab framework:
clearvars
rng(100,'twister')
uqlab

%% 2 - COMPUTATIONAL MODEL
%
% The computational model is an $8$-dimensional analytical formula 
% that is used to model the water flow through a borehole.
% The borehole function |uq_borehole| is supplied with UQLab.
%
% Create a MODEL object from the function file:
ModelOpts.mFile = 'UQLab_network_run';
myModel = uq_createModel(ModelOpts);

%%
% Type |help uq_borehole| for information on the model structure as well as
% the description of each variable.

%% 3 - PROBABILISTIC INPUT MODEL
%
% Different values for each t, w, k, n

for ii = 1:30
    InputOpts.Marginals(1,ii).Name = '$\tau$';  % time constant
    InputOpts.Marginals(1,ii).Type = 'Uniform';
    InputOpts.Marginals(1,ii).Parameters = [0.5 10];  % hr
end

for ii = 31:70
    InputOpts.Marginals(1,ii).Name = 'w';  % time constant
    InputOpts.Marginals(1,ii).Type = 'Uniform';
    InputOpts.Marginals(1,ii).Parameters = [0.8 1];  % hr
end

for ii = 71:110
    InputOpts.Marginals(1,ii).Name = 'k';  % time constant
    InputOpts.Marginals(1,ii).Type = 'Uniform';
    InputOpts.Marginals(1,ii).Parameters = [0.4 0.6];  % hr
end

for ii = 111:150
    InputOpts.Marginals(1,ii).Name = 'n';  % time constant
    InputOpts.Marginals(1,ii).Type = 'Uniform';
    InputOpts.Marginals(1,ii).Parameters = [1.4 4];  % hr
end
    
% Same values for each t, w, k, n

%    InputOpts.Marginals(1).Name = 't';  % time constant
%    InputOpts.Marginals(1).Type = 'Uniform';
%    InputOpts.Marginals(1).Parameters = [0.01 10];  % hr

%    InputOpts.Marginals(2).Name = 'w';  % reaction weight
%    InputOpts.Marginals(2).Type = 'Uniform';
%    InputOpts.Marginals(2).Parameters = [0.01 1];  % 


%    InputOpts.Marginals(3).Name = 'k';  % related to EC50 / half effect
%    InputOpts.Marginals(3).Type = 'Uniform';
%    InputOpts.Marginals(3).Parameters = [0.01 0.9];  % 

%    InputOpts.Marginals(4).Name = 'n';  % Hill coeff.
%    InputOpts.Marginals(4).Type = 'Uniform';
%    InputOpts.Marginals(4).Parameters = [1.4 5];  %


%%
% Create an INPUT object based on the specified marginals:
myInput = uq_createInput(InputOpts);

%% 4 - SENSITIVITY ANALYSIS
%
% Sobol' indices are calculated first with a direct MC simulation 
% of the model and subsequently through post-processing of the
% coefficients of its PCE and LRA approximation.

%% 4.1 MC-based Sobol' indices
%
% Select the sensitivity analysis module in UQLab
% and the Sobol' analysis method:
SobolOpts.Type = 'Sensitivity';
SobolOpts.Method = 'Sobol';

%%
% Specify the maximum order of the Sobol' indices to be calculated:
SobolOpts.Sobol.Order = 1;

%%
% Specify the sample size for the MC simulation:
SobolOpts.Sobol.SampleSize = 1000;
%%
% Note that the total cost of computation is $(M+2) \times N$,
% where $M$ is the input dimension and $N$ is the sample size.
% Therefore, the total cost for the current setup is
% $(8+2) \times 10^5 = 10^6$ evaluations of the full computational model.

%%

tic

% Run the sensitivity analysis:
mySobolAnalysisMC = uq_createAnalysis(SobolOpts);

% Check if parfor works for Sobol (https://uqworld.org/t/how-to-run-uqlab-bayesian-module-with-parallel-processing/251/3)
% parfor i = 1:152000                                                     % number of simulations
%    myAnalysisContainer{i} = uq_createAnalysis(SobolOpts,'-private'); % '-private' flag is necessary to not store the analysis object in the current UQLab session ?       
% end

toc

%%
% Retrieve the analysis results for comparison:
mySobolResultsMC = mySobolAnalysisMC.Results;


%% 5 - COMPARISON OF THE RESULTS
%
% Print the results of the Sobol' indices calculation
% based on the MC simulation:
uq_print(mySobolAnalysisMC)


%%
% Compile the relevant results for comparison:
SobolTotal = [mySobolResultsMC.Total];
SobolFirstOrder = [mySobolResultsMC.FirstOrder];

writematrix(real(SobolTotal), 'fullmodel_global_ode23s.csv') % 
writematrix(real(SobolFirstOrder), 'fullmodel_global_FO_ode23s.csv')
%% 
speciesNames = {'GLU','LPS','AGE','VEGFR1','VEGFR2','VEGF-A$_{mRNA}$','RAGE$_{ec}$','RAGE','TLR4','NADPH','NADPH$_{ec}$','ROS$_{ec}$','ROS','PI3K','AKT','PI3K$_{ec}$','AKT$_{ec}$','NF$\kappa$B$_{ec}$','NF$\kappa$B','NO','ONOO','eNOS','IL-6','TNF-$\alpha$','IL-1$\beta$','PLC-$\gamma$','VEGF-A','pJunction','Ca','Gap Width',};
reac_names = {'=\textgreater GLU','=\textgreater LPS','LPS =\textgreater TLR4','GLU =\textgreater AGE', 'AGE =\textgreater RAGE','RAGE =\textgreater NADPH','TLR4 \& ROS =\textgreater NF$\kappa$B', 'TLR4 =\textgreater PI3K','NADPH =\textgreater ROS', 'PI3K =\textgreater AKT', 'PI3K =\textgreater ROS', 'NF$\kappa$B$_{ec}$ =\textgreater TNF-$\alpha$','AKT =\textgreater NF$\kappa$B','NF$\kappa$B =\textgreater IL-6','NF$\kappa$B =\textgreater TNF-$\alpha$','NF$\kappa$B =\textgreater VEGF-A$_{mRNA}$', 'VEGF-A$_{mRNA}$ =\textgreater VEGF-A','NF$\kappa$B =\textgreater IL-1$\beta$','VEGF-A =\textgreater VEGFR1','VEGF-A =\textgreater VEGFR2','AGE =\textgreater RAGE$_{ec}$','RAGE$_{ec}$ =\textgreater NADPH$_{ec}$','VEGFR2 =\textgreater PI3K$_{ec}$','VEGFR1 =\textgreater PI3K$_{ec}$', 'NADPH$_{ec}$ =\textgreater ROS$_{ec}$', 'PI3K$_{ec}$ =\textgreater AKT$_{ec}$', 'AKT$_{ec}$ =\textgreater eNOS' , 'VEGFR1 =\textgreater PLC-$\gamma$' ,'PLC-$\gamma$ =\textgreater NF$\kappa$B$_{ec}$' , 'ROS$_{ec}$ =\textgreater NF$\kappa$B$_{ec}$','NF$\kappa$B$_{ec}$ =\textgreater IL-6', 'NF$\kappa$B$_{ec}$ =\textgreater IL-1$\beta$', 'eNOS  =\textgreater NO', 'eNOS =\textgreater ROS$_{ec}$' ,'ROS$_{ec}$ \& NO =\textgreater ONOO','!NO =\textgreater Ca','PLC-$\gamma$ =\textgreater Ca', 'Ca =\textgreater pJunction','pJunction =\textgreater Gap Width', 'Ca =\textgreater NO'};

% Create a bar plot to compare the total Sobol' indices:

% Create the plot
uq_figure('Name', 'Total Sobol'' Indices')
barWidth = 1;
figure(1)
uq_bar(1:30, real(SobolTotal(1:30,:)), barWidth)
% Set axes limits
ylim([0 1])
%xlim([0 5])
% Set labels
xlabel('Variable name')
ylabel('Total Sobol'' indices')
set(...
    gca,...
    'XTick', 1:length(InputOpts.Marginals),...
    'XTickLabel', speciesNames)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')

%% 
% Create a bar plot to compare the first-order Sobol' indices:

% Create the plot
uq_figure('Name', 'First-order Sobol'' Indices')
uq_bar(1:30, real(SobolFirstOrder), barWidth)
% Set axes limits
%xlim([0 5])
ylim([0 1])
% Set labels
xlabel('Variable name')
ylabel('First-order Sobol'' indices')
%set(...
%    gca,...
%    'XTick', 1:length(InputOpts.Marginals),...
%    'XTickLabel', [speciesNames,reac_names,reac_names]) % mySobolResultsMC.VariableNames)
% Set legend
%uq_legend({...
%    sprintf('MC-based (%.0e simulations)', mySobolResultsMC.Cost)},...
%    'Location', 'northeast')
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A', 'ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
%% 