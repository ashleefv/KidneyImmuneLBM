
speciesNames = {'GLU','LPS','AGE','VEGFR1','VEGFR2','VEGF-A$_{mRNA}$','RAGE$_{ec}$','RAGE','TLR4','NADPH','NADPH$_{ec}$','ROS$_{ec}$','ROS','PI3K','AKT','PI3K$_{ec}$','AKT$_{ec}$','NF$\kappa$B$_{ec}$','NF$\kappa$B','NO','ONOO','eNOS','IL-6','TNF-$\alpha$','IL-1$\beta$','PLC-$\gamma$','VEGF-A','pJunction','Ca','Gap Width',};
reac_names = {'=\textgreater GLU','=\textgreater LPS','LPS =\textgreater TLR4','GLU =\textgreater AGE', 'AGE =\textgreater RAGE','RAGE =\textgreater NADPH','TLR4 \& ROS =\textgreater NF$\kappa$B', 'TLR4 =\textgreater PI3K','NADPH =\textgreater ROS', 'PI3K =\textgreater AKT', 'PI3K =\textgreater ROS', 'NF$\kappa$B$_{ec}$ =\textgreater TNF-$\alpha$','AKT =\textgreater NF$\kappa$B','NF$\kappa$B =\textgreater IL-6','NF$\kappa$B =\textgreater TNF-$\alpha$','NF$\kappa$B =\textgreater VEGF-A$_{mRNA}$', 'VEGF-A$_{mRNA}$ =\textgreater VEGF-A','NF$\kappa$B =\textgreater IL-1$\beta$','VEGF-A =\textgreater VEGFR1','VEGF-A =\textgreater VEGFR2','AGE =\textgreater RAGE$_{ec}$','RAGE$_{ec}$ =\textgreater NADPH$_{ec}$','VEGFR2 =\textgreater PI3K$_{ec}$','VEGFR1 =\textgreater PI3K$_{ec}$', 'NADPH$_{ec}$ =\textgreater ROS$_{ec}$', 'PI3K$_{ec}$ =\textgreater AKT$_{ec}$', 'AKT$_{ec}$ =\textgreater eNOS' , 'VEGFR1 =\textgreater PLC-$\gamma$' ,'PLC-$\gamma$ =\textgreater NF$\kappa$B$_{ec}$' , 'ROS$_{ec}$ =\textgreater NF$\kappa$B$_{ec}$','NF$\kappa$B$_{ec}$ =\textgreater IL-6', 'NF$\kappa$B$_{ec}$ =\textgreater IL-1$\beta$', 'eNOS  =\textgreater NO', 'eNOS =\textgreater ROS$_{ec}$' ,'ROS$_{ec}$ \& NO =\textgreater ONOO','!NO =\textgreater Ca','PLC-$\gamma$ =\textgreater Ca', 'Ca =\textgreater pJunction','pJunction =\textgreater Gap Width', 'Ca =\textgreater NO'};
%%
Th_tau = real(0.1*max(max(SobolTotal([1:30],:)')));
% real(min(max(0.1*SobolTotal(1:30))));
Th_W = real(0.1*max(max(SobolTotal([31:70],:)')));
Th_k = real(0.1*max(max(SobolTotal([71:110],:)')));
Th_n = real(0.1*max(max(SobolTotal([111:150],:)')));

% Create a bar plot to compare the total Sobol' indices:

% Create the plot
uq_figure('Name', 'Total Sobol'' Indices')
barWidth = 1;
figure(1)
uq_bar(1:30, real(SobolTotal(1:30,:)), barWidth)
hold on
yline(Th_tau, ':', 'Threshold', 'LineWidth', 2)
% Set axes limits
%xlim([0 5])
% Set labels
xlabel('$\tau$')
ylabel('Total Sobol'' indices')
set(...
    gca,...
    'XTick', 1:30,...
    'XTickLabel', speciesNames)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
savefig('TotolSobol_fulltau_ode15s.fig')

figure(2)
uq_bar(31:70, real(SobolTotal(31:70,:)), barWidth)
hold on
yline(Th_W, ':', 'Threshold', 'LineWidth', 2)
% Set axes limits
%xlim([0 5])
% Set labels
xlabel('W')
ylabel('Total Sobol'' indices')
set(...
    gca,...
    'XTick', 31:70,...
    'XTickLabel', reac_names)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
savefig('TotolSobol_fullW_ode15s.fig')

figure(3)
uq_bar(71:110, real(SobolTotal(71:110,:)), barWidth)
hold on
yline(Th_k, ':', 'Threshold', 'LineWidth', 2)
% Set axes limits
%xlim([0 5])
% Set labels
xlabel('$EC_{50}$')
ylabel('Total Sobol'' indices')
set(...
    gca,...
    'XTick', 71:110,...
    'XTickLabel', reac_names)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
savefig('TotolSobol_fullK_ode15s.fig')

figure(4)
uq_bar(111:150, real(SobolTotal(111:150,:)), barWidth)
hold on
yline(Th_n, ':', 'Threshold', 'LineWidth', 2)
% Set axes limits
%xlim([0 5])
% Set labels
xlabel('n')
ylabel('Total Sobol'' indices')
set(...
    gca,...
    'XTick', 111:150,...
    'XTickLabel', reac_names)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
savefig('TotolSobol_fulln_ode15s.fig')

%%
% Create a bar plot to compare the total Sobol' indices:

% Create the plot
uq_figure('Name', 'First-order Sobol'' Indices')
barWidth = 1;
figure(5)
uq_bar(1:30, real(SobolFirstOrder(1:30,:)), barWidth)
% Set axes limits
%xlim([0 5])
% Set labels
xlabel('$\tau$')
ylabel('First-Order Sobol'' indices')
set(...
    gca,...
    'XTick', 1:30,...
    'XTickLabel', speciesNames)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
savefig("FirstOrderSobol_fulltau_ode23s.fig")

figure(6)
uq_bar(31:70, real(SobolFirstOrder(31:70,:)), barWidth)
% Set axes limits
%xlim([0 5])
% Set labels
xlabel('W')
ylabel('First Order Sobol'' indices')
set(...
    gca,...
    'XTick', 31:70,...
    'XTickLabel', reac_names)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
savefig("FirstOrderSobol_fullW_ode23s.fig")

figure(7)
uq_bar(71:110, real(SobolFirstOrder(71:110,:)), barWidth)
% Set axes limits
%xlim([0 5])
% Set labels
xlabel('$EC_{50}$')
ylabel('First-Order Sobol'' indices')
set(...
    gca,...
    'XTick', 71:110,...
    'XTickLabel', reac_names)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
savefig("FirstOrderSobol_fullk_ode23s.fig")


figure(8)
uq_bar(111:150, real(SobolFirstOrder(111:150,:)), barWidth)
% Set axes limits
%xlim([0 5])
% Set labels
xlabel('n')
ylabel('First-Order Sobol'' indices')
set(...
    gca,...
    'XTick', 111:150,...
    'XTickLabel', reac_names)
% Set legend
uq_legend({'ROS', 'IL-6', 'TNF-$\alpha$', 'IL-1$\beta$', 'VEGF-A','ROS$_{ec}$', 'NO', 'eNOS'}, 'Location', 'northeast')
savefig("FirstOrderSobol_fulln_ode23s.fig")

