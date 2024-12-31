function [Tout, Yout] = coupledODE_IVV_run(tspan, y0, params, p_params, mode, state, glu_sampled)

  global number_ctrl time_ctrl density diameter GC_conc GC_time GC_LB GC_UB time_lee glu_UB glu_LB time_finch glu_finch LB_lee UB_lee glucose_lee
  
  MC = load('data/MC_25_fen_multirun.mat');
  credible = MC.credible;

  blue = 	[0 0.4470 0.7410];


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
        for Tt = [336:1:tspan(end)]
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
    time_g = [336:1:3360];

    
    figure(3) 
    subplot(1,2,1); box;
    hold on

    plot(Tout/(24*7), Yout(:,[38]), 'LineWidth', 1.2, 'Color', 'k'); 
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
    plot(Tout/(24*7), Yout(:,37), 'LineWidth', 1.2, 'Color', 'k'); 
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
         

if state == 'diab_mice'
    
    figure(2)
    subplot(1,2,1)
    plot([336:length(GLU_p(:,1))]/(24*7), GLU_p([336:3360],1), 'LineWidth', 3, 'Color', 'k'); 
    hold on; scatter(time_finch/(7*24), glu_finch, 'MarkerFaceColor',[1 0 0],'Marker','^','Color',[1 0 0])
    hold on; errorbar(time_finch/(7*24), glu_finch, abs(glu_finch - glu_LB), abs(glu_finch - glu_UB), '^', 'Color',[1  0  0]);
    hold on; scatter(time_lee/(7*24), glucose_lee, 'MarkerFaceColor',[0  0  1],'Marker','o','Color',[0  0  1])
    hold on; errorbar(time_lee/(7*24), glucose_lee, abs(glucose_lee - LB_lee), abs(glucose_lee - UB_lee), 'o', 'Color',[0  0  1]);
    xlabel('Time (weeks)'); xlim([0,21]);
    ylabel('Glucose (mmol/l)'); ylim([0,55]);
    ax = gca; ax.FontSize = 20;
    hold on
    legend('Model', 'Finch et al. (2022)', '', 'Lee et al. (2018)', '');
    ylim([0,52]);
    x0=10;
    y0=10;
    width=800;
    height=800;
    set(gcf,'position',[x0,y0,width,height])
    ax = gca; ax.FontSize = 20;

    subplot(1,2,2)

    norm_glu(:) = (glu_finch(:) - 10) / (max(glu_UB) - 10);
    norm_UB(:) = (glu_UB(:) - 10)/(max(glu_UB) - 10);
    norm_LB(:) = (glu_LB(:) - 10)/(max(glu_UB) - 10);
    plot(Tout/(24*7), Yout(:,1), 'LineWidth', 3, 'Color', 'k');
    
    xlabel('Time (weeks)');
    ylabel('Normalized Glucose');% xlim([0,22]); 
    ax = gca; ax.FontSize = 20;
    hold on
    x0=10;
    y0=10;
    width=800;
    height=800;
    set(gcf,'position',[x0,y0,width,height])
    legend('Model')
    box on

    xlabel('Time (weeks)');
    ylabel('Normalized Glucose'); %xlim([0,22]); ylim([0,1.2]);
    ax = gca; ax.FontSize = 20;
    hold on
    legend('Normalized model')
else 
    figure(2)
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



figure(51)
for i=1:length(params{4}([1:36]))
     subplot(4,9,i)
     plot(Tout/(24*7), Yout(:,i),'LineWidth', 1.2, 'Color', 'k');
    
     hold on
     ylabel(params{4}(i))
     xlabel('Time (week)')
     xlim([0,21]);
     ax = gca; ax.FontSize = 14;
end

hold off;

end




end