function [Tout, Yout] = coupledODE_IVV_run(tspan, y0, params, p_params, mode, state, glu_sampled)

  global number_ctrl time_ctrl density diameter GC_conc GC_time GC_LB GC_UB time_lee glu_UB glu_LB time_finch glu_finch LB_lee UB_lee glucose_lee
  
  MC = load('data/MC_25_fen_multirun.mat');
  credible = MC.credible;

  blue = 	[0 0.4470 0.7410];
    start_time = 2; %weeks
    start_time_h = start_time*7*24;
    end_time = 20; %weeks
    end_time_h = end_time*7*24;

indexforNumber = 36;
indexforDiameter = 37;

  opts=[];
  %opts = odeset('RelTol',1e-20, 'MaxStep',1e-16);
  intv = "none";
  [t, y] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,params,p_params, state, glu_sampled, intv);

   Yout = real(y);
   Tout = t;

    speciesNames = params{4};

    GLU = params{1}(1,1);


    if state == 'diab_mice'
       % Gp = step_function(glu_sampled);
        for Tt = start_time_h:1:tspan(end)
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
if state == 'diab_mice'

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
    figure(40) % single run, not for publication
    % this figure runs at the output where glu_sampled = GC_conc;
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
    plot(Tout/(24*7), Yout(:,indexforNumber), 'LineWidth', 1.2, 'Color', 'k'); 
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
    plot(Tout/(24*7), Yout(:,indexforDiameter), 'LineWidth', 1.2, 'Color', 'k'); 
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

         

    Gp0 = GLU_p(start_time_h,1);
    W_GLU = (GLU_p(:,1) - Gp0) / (max(glu_UB) - Gp0);  % (1) use of max of finch et al. data error bars
    figure(2)
    figname = 'Fig2';
    yyaxis left
    plot((start_time_h:length(GLU_p(:,1)))/(24*7), GLU_p(start_time_h:end_time_h,1), 'LineWidth', 1.2, 'Color','k', 'LineStyle','--'); 
    %hold on; scatter(time_lee/(7*24), glucose_lee, 'MarkerFaceColor',[0  0  1],'Marker','o','MarkerEdgeColor',[0  0  1])
    glucose_data_Lee_sd = abs(glucose_lee - LB_lee);
    glucose_data_Finch_sd = abs(glu_finch - glu_UB); 
    hold on; errorbar(time_lee/(7*24), glucose_lee, glucose_data_Lee_sd, 'o', 'Color',[0  0  1],'MarkerFaceColor',[0  0  1]);
    % hold on; errorbar(time_lee/(7*24), glucose_lee, abs(glucose_lee - LB_lee), abs(glucose_lee - UB_lee), 'o', 'Color',[0  0  1],'MarkerFaceColor',[0  0  1]);
    %hold on; scatter(time_finch/(7*24), glu_finch, 'MarkerFaceColor',[1 0 0],'Marker','^','MarkerEdgeColor',[1 0 0])
    hold on; errorbar(time_finch/(7*24), glu_finch, glucose_data_Finch_sd, '^', 'Color',[1 0 0],'MarkerFaceColor',[1 0 0]);
    %hold on; errorbar(time_finch/(7*24), glu_finch, abs(glu_finch - glu_LB), abs(glu_finch - glu_UB), '^', 'Color',[1 0 0],'MarkerFaceColor',[1 0 0]);
    xlabel('Time (weeks)'); xlim([0,21]);
    ylabel('Glucose (mmol/l)'); 
    ax = gca; ax.FontSize = 8;
    ax.YAxis(1).Color = 'k'; 
    hold on
    
    leftymin = 4;
    leftymax = 55;
    ylim([leftymin,leftymax]);

    % need to convert the y axis to normalized units to get the scaling
    % perfect
    rightymin = (leftymin - Gp0) / (max(glu_UB) - Gp0);
    rightymax = (leftymax - Gp0) / (max(glu_UB) - Gp0);

    ax = gca; ax.FontSize = 8;

    % abs(glucose_lee - LB_lee), abs(glucose_lee - UB_lee)
    % 
    % abs(glu_finch - glu_LB), abs(glu_finch - glu_UB)

    yyaxis right
    %plot(Tout/(24*7), Yout(:,1), 'k-','LineWidth', 1.2);
    hold on
    plot((start_time_h:length(GLU_p(:,1)))/(24*7), W_GLU(start_time_h:end_time_h,1), 'LineWidth', 1.2, 'Color','k', 'LineStyle',':')
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

var = 1:37;
fig = figure(52);
figname = 'FigB';

yminSpeciesVector = zeros(1,37);
yminSpeciesVector(indexforNumber) = 4;
yminSpeciesVector(indexforDiameter) = 45;
ymaxSpeciesVector = ones(1,37);
ymaxSpeciesVector(1) = 2;
ymaxSpeciesVector(30) = 6e-3; % Actin_r
ymaxSpeciesVector(35) = 0.02; % MLC
ymaxSpeciesVector(indexforNumber) = 6.5;
ymaxSpeciesVector(indexforDiameter) = 80;

for i=var

     subplot(5,8,i)
     plot(Tout/(24*7), Yout(:,i),'LineWidth', 1.2, 'Color', 'k');

     hold on

     ylabel(params{4}(i))
     %xlabel('Time (weeks)')
     xlim([0,20]);
     ylim([yminSpeciesVector(i) ymaxSpeciesVector(i)]);
     ax = gca; ax.FontSize = 8;
end
% Common x-axis label
han = axes(fig, 'visible', 'off'); 
han.XLabel.Visible = 'on';
xlabel(han, 'Time (weeks)','FontName','Arial','FontSize',8);

widthInches = 9;
heightInches = 5;
run('ScriptForExportingImages.m')    

% in command window use ImageMagick to crop this pdf
% magick FigB.tiff FigB.pdf

end

end

if state=="norm_mice" && mode==1
    figure(41)
    plot(Tout/(24*7), Yout(:,1), 'LineWidth', 1.2);
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

    figure(42)
    subplot(1,2,1); plot(Tout/24/7, Yout(:,indexforNumber)); ylim([4,7]); ylabel("Fenestration Number")
    subplot(1,2,2); plot(Tout/24/7, Yout(:,indexforDiameter)); ylim([30,60]); ylabel("Fenestration Diameter (nm)")

    var = 1:37;
    figure(43);
    for i=var

         subplot(5,8,i)
         plot(Tout/(24*7), Yout(:,i),'LineWidth', 1.2, 'Color', 'k');
    
         hold on
    
         ylabel(params{4}(i))
         %xlabel('Time (weeks)')
         xlim([0,21]);
         ax = gca; ax.FontSize = 8;
    end
end

end