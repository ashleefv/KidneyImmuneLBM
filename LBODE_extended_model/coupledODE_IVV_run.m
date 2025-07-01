function [Tout, Yout] = coupledODE_IVV_run(tspan, y0, params, p_params, mode, state, glu_sampled)

  global number_ctrl time_ctrl density diameter GC_conc GC_time GC_LB GC_UB time_lee glu_UB glu_LB time_finch glu_finch LB_lee UB_lee glucose_lee
  
  MC = load('data/MC_25_fen_multirun.mat');
  credible = MC.credible;

  blue = 	[0 0.4470 0.7410];
    start_time = 2; %weeks
    start_time_h = start_time*7*24;
    end_time = 20; %weeks
    end_time_h = end_time*7*24;

  opts=[];
  %opts = odeset('RelTol',1e-20, 'MaxStep',1e-16);
  intv = "none";
  [t, y] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,params,p_params, state, GC_conc', intv);

   Yout = real(y);
   Tout = t;

    speciesNames = params{4};

    GLU = params{1}(1,1);


    if state == 'diab_mice'
       % Gp = step_function(glu_sampled);
        for Tt = [start_time_h:1:tspan(end)]
            if (Tt >= time_lee(1) && Tt <= time_lee(5))
                GLU_p(Tt,1) = 0.051*(Tt)  - 9.38;
                %%GLU_p(Tt,1) = -2.994e-05*Tt^2 + 0.08638*Tt - 18.4;
            else
    
                GLU_p(Tt,1) = step_function(Tt, glu_sampled);
                %GLU_p(Tt,1) = double(Gp([Tt]));
            end
        end
    
    end


%% In vivo  glucose and feedback ON


if mode == 1
  
    purple = [    0.4940    0.1840    0.5560];

 
    pop_diameter = load("data\dbmice_diameter_population.csv");
    pop_density = load("data\dbmice_density_population.csv");
    time_g = [start_time_h:1:end_time_h];

    
    figure(4) 
    figname = 'Fig4';
    subplot(2,2,1); box;
    hold on
    scatter([6], [47.91], 100, 's', 'filled', 'r');
    hold on 
    % scatter([6,10,15,20], [50.74, 60.19, 73.65, 74.63], 50, 'filled', 'k')
    % hold on
    errorbar([6,10,15,20], [50.74, 60.19, 73.65, 74.63], abs([50.74, 60.19, 73.65, 74.63] - [55.12,63.9,80.48,80]), abs([50.74, 60.19, 73.65, 74.63] - [55.12,63.9,80.48,80]), 'ko', 'MarkerFaceColor', 'k', 'LineWidth', 0.75);
    hold on
    scatter(pop_diameter(:,1), pop_diameter(:,2), 'o',  'k')
    hold on
    plot(Tout/(24*7), Yout(:,[38]), 'LineWidth', 1.2, 'Color', 'k'); 
    hold on
    [ph,msg] = jbfill(time_g/(24*7), credible(:,1,38)', credible(:,2,38)', blue, blue, 1, 0.2);
    ylabel('Fenestration Diameter (nm)'); xlabel('Time (weeks)')
%    legend('Model', '95% credible interval', 'Control data', 'Diabetes data', 'Diabetes population')
    ax = gca;  ax.FontSize = 8;    
%     t = title(['Finch et al., JASN (2022)']); t.FontSize = 8;

    subplot(2,2,2); box;
    
    hold on 
    scatter([6], [6.3], 100, 's', 'filled', 'r');
    hold on 
    % scatter([6,10,15,20], [5.65, 4.50, 4.08, 4.14], 50, 'filled', 'k')
    % hold on
    errorbar([6,10,15,20], [5.65, 4.50, 4.08, 4.14], abs([5.65, 4.50, 4.08, 4.14] - [6.6, 4.96,4.98,5.41]), abs([5.65, 4.50, 4.08, 4.14] - [4.98, 4.11, 3.15, 2.95]), 'ko', 'MarkerFaceColor', 'k',  'LineWidth', 0.75);
    hold on
    scatter(pop_density(:,1), pop_density(:,2), 'o', 'k')
    hold on
    plot(Tout/(24*7), Yout(:,37), 'LineWidth', 1.2, 'Color', 'k'); 
    hold on
    [ph,msg] = jbfill(time_g/(24*7), credible(:,1,37)', credible(:,2,37)', blue, blue, 1, 0.2);
    xlabel('Time (weeks)');
    ylabel('Fenestration Number');
    
    ax = gca;  ax.FontSize = 8;    
%     t = title(['Finch et al., JASN (2022)']); t.FontSize = 8;


    labelstring = {'A', 'B'};
    for v = 1:2
        subplot(2,2,v)
        hold on
        text(-0.15, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize',8)
        set(gca,'FontName','Arial','FontSize',8)
    end

        % overall legend
    legend( 'Control data', 'Diabetes data (mean)',  'Diabetes data (individual)','Model', '95% credible interval')
    h = legend('Location','southoutside', 'Orientation', 'horizontal');
    h.NumColumns = 3;
    p = [0.5 0.45 0.03 0.03];
    set(h,'Position', p,'Units', 'normalized');

widthInches = 6;
heightInches = 4.6;
run('ScriptForExportingImages.m')         

if state == 'diab_mice'
    Gp0 = GLU_p(start_time_h,1);
    W_GLU = (GLU_p(:,1) - Gp0) / (max(glu_UB) - Gp0);  % (1) use of max of finch et al. data error bars
    figure(2)
    figname = 'Fig2';
    yyaxis left
    co = orderedcolors("gem");
    color5 = co(5, :);
    plot([start_time_h:length(GLU_p(:,1))]/(24*7), GLU_p(start_time_h:end_time_h,1), 'LineWidth', 1.2, 'Color','k', 'LineStyle','--'); 
    %hold on; scatter(time_lee/(7*24), glucose_lee, 'MarkerFaceColor',[0  0  1],'Marker','o','MarkerEdgeColor',[0  0  1])
    hold on; errorbar(time_lee/(7*24), glucose_lee, abs(glucose_lee - LB_lee), abs(glucose_lee - UB_lee), 'o', 'Color',[0  0  1],'MarkerFaceColor',[0  0  1]);
    %hold on; scatter(time_finch/(7*24), glu_finch, 'MarkerFaceColor',[1 0 0],'Marker','^','MarkerEdgeColor',[1 0 0])
    hold on; errorbar(time_finch/(7*24), glu_finch, abs(glu_finch - glu_LB), abs(glu_finch - glu_UB), '^', 'Color',[1 0 0],'MarkerFaceColor',[1 0 0]);
    xlabel('Time (weeks)'); xlim([0,21]);
    ylabel('Glucose (mmol/l)'); 
    ax = gca; ax.FontSize = 8;
    ax.YAxis(1).Color = 'k'; 
    hold on
    
    leftymin = 5;
    leftymax = 55;
    ylim([leftymin,leftymax]);

    % need to convert the y axis to normalized units to get the scaling
    % perfect
    rightymin = (leftymin - Gp0) / (max(glu_UB) - Gp0);
    rightymax = (leftymax - Gp0) / (max(glu_UB) - Gp0);
    % x0=10;
    % y0=10;
    % width=800;
    % height=800;
    % set(gcf,'position',[x0,y0,width,height])
    ax = gca; ax.FontSize = 8;

    % abs(glucose_lee - LB_lee), abs(glucose_lee - UB_lee)
    % 
    % abs(glu_finch - glu_LB), abs(glu_finch - glu_UB)

    yyaxis right
    %plot(Tout/(24*7), Yout(:,1), 'k-','LineWidth', 1.2);
    hold on
    plot([start_time_h:length(GLU_p(:,1))]/(24*7), W_GLU(start_time_h:end_time_h,1), 'LineWidth', 1.2, 'Color','k', 'LineStyle',':')
    ax = gca; ax.FontSize = 8;
    ax.YAxis(2).Color = 'k'; 
    ylabel('Normalized Glucose Units');% xlim([0,22]);
    % rightymin = 0;
    % rightymax = 1.6;
    ylim([rightymin,rightymax]);
    box on
    legend('Simulated Glucose Input $G(t)$',   'Lee et al. (2018)', 'Finch et al. (2022)', "Normalized Glucose Reaction Weight $W^{'}_\mathrm{GLU}(t)$",'Location','southeast','interpreter','latex');
    % 
    % norm_glu(:) = (glu_finch(:) - 10) / (max(glu_UB) - 10);
    % norm_UB(:) = (glu_UB(:) - 10)/(max(glu_UB) - 10);
    % norm_LB(:) = (glu_LB(:) - 10)/(max(glu_UB) - 10);
    widthInches = 5.5;
    heightInches = 4.23;
    run('ScriptForExportingImages.m')    


else 
    figure(2)
    %figname = 'Fig2';
    plot(Tout(time_in)/(24*7), GLU_p(time_in,1), 'LineWidth', 2);
    hold on; scatter(time_finch/(7*24), ctrl_finch, 'MarkerFaceColor',[0 0 0],'Marker','o','Color',[0.75 0.5 0.25])
    hold on; errorbar(time_finch/(7*24), ctrl_finch, abs(ctrl_finch - ctrl_LB), abs(ctrl_finch - ctrl_UB), '*');
    xlabel('Time (weeks)'); %xlim([0,22]);
    ylabel('Glucose (mmol/l)'); %ylim([0,40])
    ax = gca; ax.FontSize = 8;
    hold on
    legend('Model', 'Data', '')
    x0=10;
    y0=10;
    width=800;
    height=800;
    set(gcf,'position',[x0,y0,width,height])
end


var = [1:29, 31:36]; % gap width node not plotted
fig = figure(51);
figname = 'FigB';
for i=var
    if i<=29
     subplot(4,9,i)
     plot(Tout/(24*7), Yout(:,i),'LineWidth', 1.2, 'Color', 'k');

     hold on
    end
    if i>=31
     subplot(4,9,i-1)
     plot(Tout/(24*7), Yout(:,i),'LineWidth', 1.2, 'Color', 'k');

     hold on
    end

     ylabel(params{4}(i))
     %xlabel('Time (weeks)')
     xlim([0,21]);
     ax = gca; ax.FontSize = 8;
end
% Common x-axis label
han = axes(fig, 'visible', 'off'); 
han.XLabel.Visible = 'on';
xlabel(han, 'Time (weeks)','FontName','Arial','FontSize',8);

widthInches = 9.5;
heightInches = 5;
run('ScriptForExportingImages.m')    

% in command window use ImageMagick to crop this pdf
% magick FigB.tiff FigB.pdf

end




end