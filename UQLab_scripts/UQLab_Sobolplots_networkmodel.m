function UQLab_Sobolplots_networkmodel()

%% 1 - INITIALIZE UQLAB
%
% Clear variables from the workspace,
% set random number generator for reproducible results,
% and initialize the UQLab framework:

% Must be in ~ UQLab/core directory to initialize
rng(100,'twister')
uqlab
%%
speciesNames = {'GLU','LPS','AGE','VEGFR1','VEGFR2','VEGF-A$_{mRNA}$','RAGE$_{ec}$','RAGE','TLR4','NADPH','NADPH$_{ec}$','ROS$_{ec}$','ROS','PI3K','AKT','PI3K$_{ec}$','AKT$_{ec}$','NF$\kappa$B$_{ec}$','NF$\kappa$B','NO','ONOO','eNOS','IL-6','TNF-$\alpha$','IL-1$\beta$','PLC-$\gamma$','VEGF-A','pJunction','Ca','Gap Width',};
reac_names = {'$\Rightarrow$ GLU','$\Rightarrow$ LPS','LPS $\Rightarrow$ TLR4','GLU $\Rightarrow$ AGE', 'AGE $\Rightarrow$ RAGE','RAGE $\Rightarrow$ NADPH','TLR4 \& ROS $\Rightarrow$ NF$\kappa$B', 'TLR4 $\Rightarrow$ PI3K','NADPH $\Rightarrow$ ROS', 'PI3K $\Rightarrow$ AKT', 'PI3K $\Rightarrow$ ROS', 'NF$\kappa$B$_{ec}$ $\Rightarrow$ TNF-$\alpha$','AKT $\Rightarrow$ NF$\kappa$B','NF$\kappa$B $\Rightarrow$ IL-6','NF$\kappa$B $\Rightarrow$ TNF-$\alpha$','NF$\kappa$B $\Rightarrow$ VEGF-A$_{mRNA}$', 'VEGF-A$_{mRNA}$ $\Rightarrow$ VEGF-A','NF$\kappa$B $\Rightarrow$ IL-1$\beta$','VEGF-A $\Rightarrow$ VEGFR1','VEGF-A $\Rightarrow$ VEGFR2','AGE $\Rightarrow$ RAGE$_{ec}$','RAGE$_{ec}$ $\Rightarrow$ NADPH$_{ec}$','VEGFR2 $\Rightarrow$ PI3K$_{ec}$','VEGFR1 $\Rightarrow$ PI3K$_{ec}$', 'NADPH$_{ec}$ $\Rightarrow$ ROS$_{ec}$', 'PI3K$_{ec}$ $\Rightarrow$ AKT$_{ec}$', 'AKT$_{ec}$ $\Rightarrow$ eNOS' , 'VEGFR1 $\Rightarrow$ PLC-$\gamma$' ,'PLC-$\gamma$ $\Rightarrow$ NF$\kappa$B$_{ec}$' , 'ROS$_{ec}$ $\Rightarrow$ NF$\kappa$B$_{ec}$','NF$\kappa$B$_{ec}$ $\Rightarrow$ IL-6', 'NF$\kappa$B$_{ec}$ $\Rightarrow$ IL-1$\beta$', 'eNOS  $\Rightarrow$ NO', 'eNOS $\Rightarrow$ ROS$_{ec}$' ,'ROS$_{ec}$ \& NO $\Rightarrow$ ONOO','!NO $\Rightarrow$ Ca','PLC-$\gamma$ $\Rightarrow$ Ca', 'Ca $\Rightarrow$ pJunction','pJunction $\Rightarrow$ Gap Width', 'Ca $\Rightarrow$ NO'};

%% read stored Sobol coefficients in data files

SobolTotal = readmatrix('fullmodel_global_Total__ode23s.csv'); % 
SobolFirstOrder = readmatrix('fullmodel_global_FirstOrder_ode23s.csv');

%%

Th_tau = real(0.1*max(max(SobolTotal([1:30],:)')));
Th_W = real(0.1*max(max(SobolTotal([31:70],:)')));
Th_k = real(0.1*max(max(SobolTotal([71:110],:)')));
Th_n = real(0.1*max(max(SobolTotal([111:150],:)')));

% Create a bar plot to compare the total Sobol' indices:

% Create the plot
uq_figure('Name', 'Total Sobol'' Indices')
barWidth = 1;


figure(5)
uq_bar(1:30, real(SobolTotal(1:30,:)), barWidth)
hold on
yline(Th_tau, ':', 'Threshold', 'LineWidth', 2)
% Set axes limits
%xlim([0 5])
% Set labels
xlabel('Time constant ($\tau$)')
ylabel('Total Sobol'' indices')
set(...
    gca,...
    'FontSize',14, ...
    'XTick', 1:30,...
    'XTickLabel', speciesNames)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
%savefig('TotolSobol_fulltau_ode23s.fig')
grid off

figure(6)
uq_bar(31:70, real(SobolTotal(31:70,:)), barWidth)
hold on
yline(Th_W, ':', 'Threshold', 'LineWidth', 2)
% Set axes limits
%xlim([0 5])
% Set labels
xlabel('Reaction weight ($W$) for reaction rules')
ylabel('Total Sobol'' indices')
set(...
    gca,...
    'FontSize',14, ...
    'XTick', 31:70,...
    'XTickLabel', reac_names)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
%savefig('TotolSobol_fullW_ode23s.fig')
grid off

figure(7)
uq_bar(71:110, real(SobolTotal(71:110,:)), barWidth)
hold on
yline(Th_k, ':', 'Threshold', 'LineWidth', 2)
% Set axes limits
%xlim([0 5])
% Set labels
xlabel('Half effect ($EC_{50}$) for reaction rules')
ylabel('Total Sobol'' indices')
set(...
    gca,...
    'FontSize',14, ...
    'XTick', 71:110,...
    'XTickLabel', reac_names)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
%savefig('TotolSobol_fullK_ode23s.fig')
grid off

figure(8)
uq_bar(111:150, real(SobolTotal(111:150,:)), barWidth)
hold on
yline(Th_n, ':', 'Threshold', 'LineWidth', 2)
% Set axes limits
%xlim([0 5])
% Set labels
xlabel('Hill coefficient ($n$) for reaction rules')
ylabel('Total Sobol'' indices')
set(...
    gca,...
    'FontSize',14, ...
    'XTick', 111:150,...
    'XTickLabel', reac_names)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
%savefig('TotolSobol_fulln_ode23s.fig')
grid off

%%
% Create a bar plot to compare the total Sobol' indices:

% Create the plot
uq_figure('Name', 'First-order Sobol'' Indices')
barWidth = 1;
figure(1)
uq_bar(1:30, real(SobolFirstOrder(1:30,:)), barWidth)
% Set axes limits
%xlim([0 5])
% Set labels
xlabel('Time constant ($\tau$)')
ylabel('First-Order Sobol'' indices')
set(...
    gca,...
    'FontSize',14, ...
    'XTick', 1:30,...
    'XTickLabel', speciesNames)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
%savefig("FirstOrderSobol_fulltau_ode23s.fig")
grid off

figure(2)
uq_bar(31:70, real(SobolFirstOrder(31:70,:)), barWidth)
% Set axes limits
%xlim([0 5])
% Set labels
xlabel('Reaction weight ($W$) for reaction rules')
ylabel('First Order Sobol'' indices')
set(...
    gca,...
    'FontSize',14, ...
    'XTick', 31:70,...
    'XTickLabel', reac_names)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
%savefig("FirstOrderSobol_fullW_ode23s.fig")
grid off

figure(3)
uq_bar(71:110, real(SobolFirstOrder(71:110,:)), barWidth)
% Set axes limits
%xlim([0 5])
% Set labels
xlabel('Half effect ($EC_{50}$) for reaction rules')
ylabel('First-Order Sobol'' indices')
set(...
    gca,...
    'FontSize',14, ...
    'XTick', 71:110,...
    'XTickLabel', reac_names)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
%savefig("FirstOrderSobol_fullk_ode23s.fig")
grid off

figure(4)
uq_bar(111:150, real(SobolFirstOrder(111:150,:)), barWidth)
% Set axes limits
%xlim([0 5])
% Set labels
xlabel('Hill coefficient ($n$) for reaction rules')
ylabel('First-Order Sobol'' indices')
set(...
    gca,...
    'FontSize',14, ...
    'XTick', 111:150,...
    'XTickLabel', reac_names)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
%savefig("FirstOrderSobol_fulln_ode23s.fig")
grid off

end
