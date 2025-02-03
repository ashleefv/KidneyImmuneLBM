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
speciesNames = {'1: GLU','2: LPS','3: AGE','4: VEGFR1','5: VEGFR2','6: VEGF-A$_{mRNA}$','7: RAGE$_{ec}$','8: RAGE','9: TLR4','10: NADPH','11: NADPH$_{ec}$','12: ROS$_{ec}$','13: ROS','14: PI3K','15: AKT','16: PI3K$_{ec}$','17: AKT$_{ec}$','18: NF$\kappa$B$_{ec}$','19: NF$\kappa$B','20: NO','21: ONOO','22: eNOS','23: IL-6','24: TNF-$\alpha$','25: IL-1$\beta$','26: PLC-$\gamma$','27: VEGF-A','28: pJunction','29: Ca',};
reac_names = {'1: $\Rightarrow$ GLU','2: $\Rightarrow$ LPS','3: LPS $\Rightarrow$ TLR4','4: GLU $\Rightarrow$ AGE', '5: AGE $\Rightarrow$ RAGE','6: RAGE $\Rightarrow$ NADPH','7: TLR4 \& ROS $\Rightarrow$ NF$\kappa$B', '8: TLR4 $\Rightarrow$ PI3K','9: NADPH $\Rightarrow$ ROS', '10: PI3K $\Rightarrow$ AKT', '11: PI3K $\Rightarrow$ ROS', '12: NF$\kappa$B$_{ec}$ $\Rightarrow$ TNF-$\alpha$','13: AKT $\Rightarrow$ NF$\kappa$B','14: NF$\kappa$B $\Rightarrow$ IL-6','15: NF$\kappa$B $\Rightarrow$ TNF-$\alpha$','16: NF$\kappa$B $\Rightarrow$ VEGF-A$_{mRNA}$', '17: VEGF-A$_{mRNA}$ $\Rightarrow$ VEGF-A','18: NF$\kappa$B $\Rightarrow$ IL-1$\beta$','19: VEGF-A $\Rightarrow$ VEGFR1','20: VEGF-A $\Rightarrow$ VEGFR2','21: AGE $\Rightarrow$ RAGE$_{ec}$','22: RAGE$_{ec}$ $\Rightarrow$ NADPH$_{ec}$','23: VEGFR2 $\Rightarrow$ PI3K$_{ec}$','24: VEGFR1 $\Rightarrow$ PI3K$_{ec}$', '25: NADPH$_{ec}$ $\Rightarrow$ ROS$_{ec}$', '26: PI3K$_{ec}$ $\Rightarrow$ AKT$_{ec}$', '27: AKT$_{ec}$ $\Rightarrow$ eNOS' , '28: VEGFR1 $\Rightarrow$ PLC-$\gamma$' ,'29: PLC-$\gamma$ $\Rightarrow$ NF$\kappa$B$_{ec}$' , '30: ROS$_{ec}$ $\Rightarrow$ NF$\kappa$B$_{ec}$','31: NF$\kappa$B$_{ec}$ $\Rightarrow$ IL-6', '32: NF$\kappa$B$_{ec}$ $\Rightarrow$ IL-1$\beta$', '33: eNOS  $\Rightarrow$ NO', '34: eNOS $\Rightarrow$ ROS$_{ec}$' ,'35: ROS$_{ec}$ \& NO $\Rightarrow$ ONOO','36: !NO $\Rightarrow$ Ca','37: PLC-$\gamma$ $\Rightarrow$ Ca', '38: Ca $\Rightarrow$ pJunction','39: Ca $\Rightarrow$ NO'};

%% read stored Sobol coefficients in data files

SobolTotal = readmatrix('fullmodel_global_Total__ode23s.csv'); % 
SobolFirstOrder = readmatrix('fullmodel_global_FirstOrder_ode23s.csv');

%%

Th_tau = real(0.1*max(max(SobolTotal([1:29],:)')));
Th_W = real(0.1*max(max(SobolTotal([31:70],:)')));
Th_k = real(0.1*max(max(SobolTotal([71:110],:)')));
Th_n = real(0.1*max(max(SobolTotal([111:150],:)')));

% Create a bar plot to compare the total Sobol' indices:

% Create the plot
uq_figure('Name', 'Total Sobol'' Indices')
barWidth = 1;


figure(5)
uq_bar(1:29, real(SobolTotal(1:29,:)), barWidth)
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
    'XTick', 1:29,...
    'XTickLabel', speciesNames)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
%savefig('TotolSobol_fulltau_ode23s.fig')
grid off

figure(6)
uq_bar([31:68,70], real(SobolTotal([31:68,70],:)), barWidth)
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
    'XTick', [31:69],...
    'XTickLabel', reac_names)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
%savefig('TotolSobol_fullW_ode23s.fig')
grid off

figure(7)
uq_bar([71:108,110], real(SobolTotal([71:108,110],:)), barWidth)
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
    'XTick', [71:109],...
    'XTickLabel', reac_names)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
%savefig('TotolSobol_fullK_ode23s.fig')
grid off

figure(8)
uq_bar([111:148,150], real(SobolTotal([111:148,150],:)), barWidth)
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
    'XTick', [111:149],...
    'XTickLabel', reac_names)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
%savefig('TotolSobol_fulln_ode23s.fig')
grid off

%%
% Create a bar plot to compare the total Sobol' indices:
SobolFirstOrder2 = SobolFirstOrder; % duplicate
SobolFirstOrder2([30,69,109,149],:) = []; % columns removed
% Create the plot
uq_figure('Name', 'First-order Sobol'' Indices')
barWidth = 1;
figure(1)
uq_bar(1:29, real(SobolFirstOrder2(1:29,:)), barWidth)
hold on
% Set axes limits
%xlim([0 5])
% Set labels
xlabel('Time constant ($\tau$)')
ylabel('First-Order Sobol'' indices')
set(...
    gca,...
    'FontSize',14, ...
    'XTick', 1:29,...
    'XTickLabel', speciesNames)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
%savefig("FirstOrderSobol_fulltau_ode23s.fig")
grid off

figure(2)
uq_bar([30:68], real(SobolFirstOrder2([30:68],:)), barWidth)
% Set axes limits
%xlim([0 5])
% Set labels
xlabel('Reaction weight ($W$) for reaction rules')
ylabel('First Order Sobol'' indices')
set(...
    gca,...
    'FontSize',14, ...
    'XTick', [30:68],...
    'XTickLabel', reac_names)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
%savefig("FirstOrderSobol_fullW_ode23s.fig")
grid off

figure(3)
uq_bar([69:107], real(SobolFirstOrder2([69:107],:)), barWidth)
% Set axes limits
%xlim([0 5])
% Set labels
xlabel('Half effect ($EC_{50}$) for reaction rules')
ylabel('First-Order Sobol'' indices')
set(...
    gca,...
    'FontSize',14, ...
    'XTick', [69:107],...
    'XTickLabel', reac_names)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
%savefig("FirstOrderSobol_fullk_ode23s.fig")
grid off

figure(4)
uq_bar([108:146], real(SobolFirstOrder2([108:146],:)), barWidth)
% Set axes limits
%xlim([0 5])
% Set labels
xlabel('Hill coefficient ($n$) for reaction rules')
ylabel('First-Order Sobol'' indices')
set(...
    gca,...
    'FontSize',14, ...
    'XTick', [108:146],...
    'XTickLabel', reac_names)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
%savefig("FirstOrderSobol_fulln_ode23s.fig")
grid off

end
