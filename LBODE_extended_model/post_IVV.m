function [s_FD_Ym, s_FD_W] = post_IVV(params, y0, tspan, p_params, state, task, Tstop)

global number_ctrl time_ctrl density diameter GC_conc GC_time GC_LB GC_UB time_lee glu_UB glu_LB time_finch glu_finch LB_lee UB_lee glucose_lee

intv = 'none';
glu_sampled = zeros(11,1);

for i = 1:length(GC_time)
           glu_sampled(i) = unifrnd(GC_LB(:,i), GC_UB(:,i)); % 
end


RP = length(params{1}(1,:));
SP = length(params{3}(:));

s_FD_Ym = []; s_FD_W = [];
%% 1: test_knockout
if task == 1
    opts=[];
    Nn = 10;
    state = 'diab_mice';
    
    glu_sampled = zeros(11,Nn);

    for Nstep = 1:Nn
        for i = 1:length(GC_time)
            glu_sampled(i,Nstep) = unifrnd(GC_LB(:,i), GC_UB(:,i)); %

        end

        [t, y] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,params,p_params, state, glu_sampled(:,Nstep), intv);
        YstepP(Nstep,:) = real(y(end,:));
        disp(size(YstepP))
    end

    s_ctrl_n([1:Nn],1,1) = 5.8; % reference is set for initial data value for healthy case
    s_ctrl_d([1:Nn],1,1) = 50;  % reference is set for initial data value for healthy case
    test_i = [29,33,32,35,2]; % inhibitors: KN93, ML7, Y27632, CalA, CytB
    z_params = params;
    p_n = zeros(1,6); p_d = zeros(1,6); p_c_n = zeros(1,7); p_c_d = zeros(1,7);

    for inh = [1:length(test_i)]
        z_params{3}(test_i(inh)) = 0;
        

        for Nstep = 1:Nn
            
            [tn, yn] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,z_params,p_params, state, glu_sampled(:,Nstep), intv);
            YstepP_new(Nstep,:) = real(yn(end,:));
            fprintf('run %i finished\n', Nstep)

        end
        
        disp(size(YstepP_new))

    
%        figure(29); hold on;
%        subplot(2,1,1); hold on; plot(tn/(24*7), yn(:,37), 'LineWidth', 2, 'LineStyle', linest(inh)); xlabel('Time (week)'); ylabel('Fenestration Number')
%        subplot(2,1,2); hold on; plot(tn/(24*7), yn(:,38), 'LineWidth', 2, 'LineStyle', linest(inh)); xlabel('Time (week)'); ylabel('Fenestration Diameter') 
        
        s_FC(:,1,:) = YstepP(:,:); % decide what should be control
       
        s_FC(:,inh+1,:) = YstepP_new(:,:);

   
       
        [h, p_n(1,inh+1)] = ttest2(s_FC(:,inh+1,37), s_FC(:,1,37), 'Alpha', 0.05);
        [h, p_d(1,inh+1)] = ttest2(s_FC(:,inh+1,38), s_FC(:,1,38), 'Alpha', 0.05);
        [h, p_c_n(1,inh)] = ttest2(s_FC(:,inh,37), s_ctrl_n(1,1,:), 'Alpha', 0.05);
        [h, p_c_d(1,inh)] = ttest2(s_FC(:,inh,38), s_ctrl_d(1,1,:), 'Alpha', 0.05);
        z_params = params;
        
    end
    

    s_FC_d = squeeze(s_FC(:,:,38));
    s_FC_n = squeeze(s_FC(:,:,37));

    CI_diameter_ub = mean(s_FC_d,1) + std(s_FC_d,1)*1.96; % upper-bound for diameter predicted
    CI_diameter_lb = mean(s_FC_d,1) - std(s_FC_d,1)*1.96;

    CI_number_ub = mean(s_FC_n,1) + std(s_FC_n,1)*1.96; % upper-bound for number predicted
    CI_number_lb = mean(s_FC_n,1) - std(s_FC_n,1)*1.96; 

    % disp(mean(s_FC_d,1)); disp(std(s_FC_d,1));
    % disp(mean(s_FC_n,1)); disp(std(s_FC_n,1));

%%
     
     figure(6); subplot(2,1,1); b = bar([mean(s_ctrl_n), mean(s_FC([1:Nn],:,37))], 'white'); xticks([1:length(test_i)+2]); xticklabels({'healthy', 'diseased (no treatment)', 'KN93', 'ML7', 'Y27632', 'CalA', 'CytB'}); ylabel('Fenestration Number');
     b.FaceColor = 'flat';
     b.CData(1,:) = [0 0 1];
     b.CData(2,:) = [0 0 0];
     b.CData(3,:) = [1 1 1]; b.CData(4,:) = [1 1 1]; b.CData(5,:) = [1 1 1]; b.CData(6,:) = [1 1 1]; b.CData(7,:) = [1 1 1];
   hold on; subplot(2,1,1); er = errorbar([2:7], [mean(s_FC([1:Nn],:,37))], (CI_number_lb-CI_number_ub)); er.Color = 'r';          er.LineStyle = 'none'; er.LineWidth=2;  


     
     groups1 = {[2,5],[2,6],[2,7]}; H=sigstar(groups1,[0,0,0]);
     groupsc = {[1,2], [1,3], [1,4], [1,5], [1,6], [1,7]}; H=sigstar(groupsc,[0,0,0,0,0.0002,0.189]*1e-14); set(H,'color','b')
%    subplot(2,1,1); legend('','*** p_{value}<1E-3', '', '', '*** p_{value}<1E-3')
     set(gca,'FontSize',12)

     figure(6); subplot(2,1,2); B = bar([mean(s_ctrl_d), mean(s_FC([1:Nn],:,38))], 'white'); xticks([1:length(test_i)+2]); xticklabels({'healthy', 'diseased (no treatment)', 'KN93', 'ML7', 'Y27632', 'CalA', 'CytB'}); ylabel('Fenestration Diameter (nm)')
     B.FaceColor = 'flat';
     B.CData(1,:) = [0 0 1];
     B.CData(2,:) = [0 0 0];
     B.CData(3,:) = [1 1 1]; B.CData(4,:) = [1 1 1]; B.CData(5,:) = [1 1 1]; B.CData(6,:) = [1 1 1]; B.CData(7,:) = [1 1 1];

     hold on; subplot(2,1,2); er = errorbar([2:7], [mean(s_FC([1:Nn],:,38))], (CI_diameter_lb - CI_diameter_ub)); er.Color = 'r';  er.LineStyle = 'none'; er.LineWidth=2; 
     groups2 = {[2,5]}; H=sigstar(groups2,[0]);
     groupsC = {[1,2], [1,3], [1,4], [1,5], [1,6], [1,7]}; H=sigstar(groupsC,[0,0.0061, 0.1679, 0, 0.0662, 0.0007]*1e-11); set(H,'color','b')

     set(gca,'FontSize',12)


     %%
     % figure(6); subplot(2,1,1); b = bar([mean(s_FC([1:Nn],:,37))], 'white'); xticks([1:length(test_i)+2]); xticklabels({'healthy', 'KN93', 'ML7', 'Y27632', 'CalA', 'CytB'}); ylabel('Fenestration Number');
     % b.FaceColor = 'flat';
     % b.CData(1,:) = [0 0 1];
     % b.CData(2,:) = [1 1 1];
     % b.CData(3,:) = [1 1 1]; b.CData(4,:) = [1 1 1]; b.CData(5,:) = [1 1 1]; b.CData(6,:) = [1 1 1]; 
     % 
     % set(gca,'FontSize',12)
     % 
     % figure(6); subplot(2,1,2); B = bar([mean(s_FC([1:Nn],:,38))], 'white'); xticks([1:length(test_i)+2]); xticklabels({'healthy',  'KN93', 'ML7', 'Y27632', 'CalA', 'CytB'}); ylabel('Fenestration Diameter (nm)')
     % B.FaceColor = 'flat';
     % B.CData(1,:) = [0 0 1];
     % B.CData(2,:) = [1 1 1];
     % B.CData(3,:) = [1 1 1]; B.CData(4,:) = [1 1 1]; B.CData(5,:) = [1 1 1]; B.CData(6,:) = [1 1 1]; 
     % 
     % set(gca,'FontSize',12)


end
%% 2: LSA-based perturbation
if task ==2

% sensitivity coefficient initialization
% size: length(species), length(params), length(time)
s_FD_Ym = zeros(SP, SP, 1); 
s_FD_W = zeros(SP, RP, 1);

linest = ["--", "-.", ":", "--", "-.", ":",  "--", "-.", ":", "--", "-.", ":", "--", "-.", ":"];

tspan = [336:1:3360];
percent = 50; opts=[];


[t, y] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,params,p_params, state, GC_conc', intv); 
y = abs(y);

% Parameter: Ymax
%deltaP = -1; % full knockdown

for m = 1:SP
    params_new = params;
    params_new{3}(m) = 0; %params{3}(m)*(1 + deltaP); % perturb each "Ymax" parameter by a small amount
    opts = [];
    [time,dy_model] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,params_new,p_params, state, GC_conc', intv);
    

    dy_model = real(dy_model);

    dy_modelR = real(dy_model);
    
    
    for l = 1:size(dy_modelR,2)
        s_FD_Ym(l,m,:) = (dy_modelR(end,l) - y(end,l))/(params_new{3}(m) - params{3}(m)); %/(percent*1e-2); % at 20 weeks
    end
%     figure(2);
%     subplot(5,7,m);
%     plot(time/(24*7), dy_modelR(:,37)); hold on; legend(params{4}(m))
end
% gap width removed from plotting (node 30)
figure(53); subplot(2,1,1); heatmap(real(s_FD_Ym([37],[1:29,31:36],1)), 'Colormap', jet); ax = gca; ax.XData = params{4}([1:29,31:36]); ax.YData = params{4}([37])  ; ax.FontSize = 12; 
figure(54); subplot(2,1,1); heatmap(real(s_FD_Ym([38],[1:29,31:36],1)), 'Colormap', jet); ax = gca; ax.XData = params{4}([1:29,31:36]); ax.YData = params{4}([38])  ; ax.FontSize = 12; 


%%
% Parameter: W
for m = 1:RP
    params_new = params;
    params_new{1}(1,m) = 0; %params{1}(1,m)*(1 + deltaP); % perturb each "W" parameter by a small amount
        
    [time,dy_model] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,params_new,p_params, state, GC_conc', intv);
     
    dy_model = real(dy_model);
    
    
    dy_modelR = real(dy_model);
    
    for l = 1:size(dy_modelR,2)
        s_FD_W(l,m,:) = ((dy_modelR(end,l)) - y(end,l))/(params_new{1}(1,m) - params{1}(1,m)); %/(percent*1e-2); % difference in values at 20 weeks
    end
end

figure(53); subplot(2,1,2); heatmap(real(s_FD_W([37],[2:29,31:36],1)), 'Colormap', jet); ax = gca; ax.XData = params{5}([2:29,31:36]);  ax.YData = params{4}([37]);   ax.FontSize = 12; 

figure(54); subplot(2,1,2); heatmap(real(s_FD_W([38],[2:29,31:36],1)), 'Colormap', jet); ax = gca; ax.XData = params{5}([2:29,31:36]);  ax.YData = params{4}([38]);  ax.FontSize = 12; 
% figure(19); bar(real(s_FD_W(37,:,1)));
% set(gca, 'XTick',1:length(RP), 'XTickLabel',params{5}(:))
% grid on;

% z_params = params;
% Pin = [0:0.1:1];
% inhib_prod = [2,4,17,20,21,22,25,30,32,43,45,48]; %[4,16,17,20,21,26,27,28,37,41,42,43,45]; 3,16 removed
% inhib_knock =  [1,3,5,6,7,11,12,18,19,27,32,34]; % [2,31,33,34,35,36]; %9,25 removed
% 
%     for inh=[1:length(inhib_knock)]
%         z_params = params;
%         
%         for pinh = [1:length(Pin)]
%         
%             z_params{3}(inhib_knock(inh)) = Pin(pinh);
%         
%             opts=[];
%             [T, Y] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,z_params,p_params, state, glu_sampled, intv);
%             Y = real(Y);
%             
%         
%             Y_fc(inh,pinh) = real(Y(end,37));
%             Y_d(inh,pinh)  = real(Y(end,38));
%     
%             
%         end
%         perc(inh,1) = (1 - Pin(find(Y_d(inh,:)<50,1,'last')))*100;
%         perc(inh,2) = (1 - Pin(find(Y_fc(inh,:)>6,1,'last')))*100;
% 
%         
%         %subplot(3,5,inh); scatter(Pin, Y_d(inh,:), '*'); xlabel(append('y_{max}(', params{4}(inh), ')')); ylabel('Diameter')
%     end
%     figure(21)
%     bar(perc(:,:)); xticks([1:length(inhib_knock)]); xticklabels(params{4}(inhib_knock))
% 
% 
% 
% 
% 
% z_params = params;
% 
% for inh=[1:length(inhib_knock)]
%     z_params = params;
%     for pinh = [1:length(Pin)]
%     
%     %z_params{3}(inhib_knock(inh)) = 0;
%     z_params{3}(inhib_knock(inh)) = Pin(pinh);
% 
%     opts=[];
%     [T, Y] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,z_params,p_params, state, glu_sampled, intv);
%     Y = real(Y);
%     FenC(1) =  y0(37);
%     for Tt = [2:1:length(T)]
%         FenC(Tt) = (1/(1 + Y(Tt,2))) + (Y(Tt,31)/(1 + Y(Tt,31))) - (1/(1 + Y(Tt,2)))*(Y(Tt,31)/(1 + Y(Tt,31)));
%     
%     end
%     Y(:,37) = FenC(:);
% 
%     %disp(real(Y(end,37)) - real(y(end,37)));
%     %z_params{3}(inh) = params{3}(inh);
%  
%     Y_fc(inh,pinh) = real(Y(end,37));
%     Y_d(inh,pinh)  = real(Y(end,38));
% 
% 
% %     figure(28); hold on;
% %     subplot(2,2,1); hold on; plot(T/(24*7), real(Y(:,37)), linest(inh), 'LineWidth', 2); 
% %     subplot(2,2,2); hold on; plot(T/(24*7), real(Y(:,38)), linest(inh), 'LineWidth', 2); 
% 
%     
%     end
%     
% %     figure(21)
% %     subplot(3,5,inh); scatter(Pin, Y_d(inh,:), '*'); xlabel(append('y_{max}(', params{4}(inh), ')')); ylabel('Diameter')
% end
%
% figure(28); hold on; subplot(2,2,1); lgd = legend(['default', params{4}(inhib_knock)], 'Location', 'NorthEast'); xlabel('Time (weeks)'); ylabel('Normalized Change in Fenestration Number');
% figure(28); hold on; subplot(2,2,2); lgd = legend(['default', params{4}(inhib_knock)], 'Location', 'NorthEast'); xlabel('Time (weeks)'); ylabel('Normalized Change in Fenestration Diameter');
% lgd.FontSize = 14;
% 
% z_params = params;
% for inh=[1:length(inhib_prod)]
%     z_params = params;
% %     for pinh = [1:length(Pin)]
%     
%     z_params{1}(1, inhib_prod(inh)) = 0;
% %     z_params{1}(1,inhib_prod(inh)) = Pin(pinh);
% 
%     opts=[];
%     [T, Y] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,z_params,p_params, state, glu_sampled, intv);
%     Y = real(Y);
%     FenC(1) =  y0(37); 
%     A = 2.06;
%     for Tt = [2:1:length(T)]
%         FenC(Tt) = (1/(1 + A*Y(Tt,2))) + (Y(Tt,31)/(1 + Y(Tt,31))) - 2*(1/(1 + A*Y(Tt,2)))*(Y(Tt,31)/(1 + Y(Tt,31)));
%     end
%     Y(:,37) = FenC(:);
% 
% %     Y_fc(inh,pinh) = real(Y(end,37));
% %     Y_d(inh,pinh)  = real(Y(end,38));
%     figure(28); hold on;
%     subplot(2,2,3); hold on; plot(T/(24*7), real(Y(:,37)), linest(inh), 'LineWidth', 2); 
%     subplot(2,2,4); hold on; plot(T/(24*7), real(Y(:,38)), linest(inh), 'LineWidth', 2); 
%     %z_params{1}(1, inh) = params{1}(1, inh);
%   
% %     end
% %     figure(21)
% %     subplot(3,5,inh); scatter(Pin, Y_d(inh,:), '*'); xlabel(append('W(', params{5}(inh), ')')); ylabel('Diameter')
% end
% figure(28); hold on; subplot(2,2,3); lgd = legend(['default', params{5}(inhib_prod)], 'Location', 'NorthEast'); xlabel('Time (weeks)'); ylabel('Normalized Change in Fenestration Number');
% figure(28); hold on; subplot(2,2,4); lgd = legend(['default', params{5}(inhib_prod)], 'Location', 'NorthEast'); xlabel('Time (weeks)'); ylabel('Normalized Change in Fenestration Diameter');
% lgd.FontSize = 14;


end

%%  3: time-dependent intervention
if task == 3

    % sensitivity coefficient initialization
    % size: length(species), length(params), length(time)
    s_FD_Ym = zeros(SP, SP, 1); 
    s_FD_W = zeros(SP, RP, 1);
    
    linest = ["--", "-.", ":", "--", "-.", ":",  "--", "-.", ":", "--", "-.", ":", "--", "-.", ":"];
    
    tspan = [336:1:Tstop];
    percent = 50; opts=[];
    
    [t, y] = ode15s(@coupledODE_IVV_step,tspan,y0,opts,params,p_params, state, GC_conc', intv); 
    y = abs(y);

    [tfull, yfull] = ode15s(@coupledODE_IVV_step,[336:5000],y0,opts,params,p_params, state, GC_conc', intv); 
    yfull = abs(yfull);
    
    
    
    figure(Tstop); hold on; ax = gca; ax.FontSize = 14;
    subplot(2,2,1); box; hold on; plot(tfull/(24*7), yfull(:,37), 'LineWidth', 2, 'color', 'k'); 
    subplot(2,2,2); box; hold on; plot(tfull/(24*7), yfull(:,37), 'LineWidth', 2, 'color', 'k'); 
    subplot(2,2,3); box; hold on; plot(tfull/(24*7), yfull(:,38), 'LineWidth', 2, 'color', 'k'); 
    subplot(2,2,4); box; hold on; plot(tfull/(24*7), yfull(:,38), 'LineWidth', 2, 'color', 'k'); 
    
    
    inhib_knock_n = [1,2,3,5,6,7,9,11,12,18,19,25,27,32,34]; % [2,31,33,34,35,36]; %9,25 removed
    inhib_prod_n = [47,21,24,16,17,2,3,18,20,29,42,44];
    % % inhib_knock_n2 = [35,36,31]; %promotes fenestrations
    
    tnew = [Tstop:1:5000];
    y0(:) = y(end,:);
    z_params = params;
    for inh = [1:length(inhib_knock_n)]
        z_params{3}(inhib_knock_n(inh)) = 0.5;
        [tn, yn] = ode15s(@coupledODE_IVV_step,tnew,y0,opts,z_params,p_params, state, GC_conc', intv); 
        yn = abs(yn);
        
        
        figure(Tstop); hold on;
        subplot(2,2,1); hold on; plot(tn/(24*7), yn(:,37), 'LineWidth', 2, 'LineStyle', linest(inh)); xlabel('Time (week)'); ylabel('Fenestration Number')
        
        z_params = params;
    end

    lgd = legend(['default', params{4}(inhib_knock_n)], 'Location', 'NorthEast');

    z_params = params;
    for inh = [1:length(inhib_prod_n)]
        %z_params{3}(inhib_knock_n2(inh)) = 2;
        z_params{1}(1,inhib_prod_n(inh)) = 0.5;
        [tn, yn] = ode15s(@coupledODE_IVV_step,tnew,y0,opts,z_params,p_params, state, GC_conc', intv); 
        yn = abs(yn);
        
        
        figure(Tstop); hold on; ax = gca; ax.FontSize = 14;
        
        subplot(2,2,2); hold on; plot(tn/(24*7), yn(:,37), 'LineWidth', 2, 'LineStyle', linest(inh)); xlabel('Time (week)'); ylabel('Fenestration Number')
        
        z_params = params;
    end
    lgd = legend(['default', params{5}(inhib_prod_n)], 'Location', 'NorthEast');
    
    inhib_knock_d = [1,3,5,6,7,9,11,12,18,19,25,27,32,34];
    inhib_prod = [2,3,19,20,29,42,44,15,24,21];
   
    
    z_params = params;
    for inh = [1:length(inhib_prod)]
    z_params{1}(1,inhib_prod(inh)) = 0.5;
    [tn, yn] = ode15s(@coupledODE_IVV_step,tnew,y0,opts,z_params,p_params, state, GC_conc', intv); 
    yn = abs(yn);
    
    
    figure(Tstop); hold on;
    subplot(2,2,4); hold on; plot(tn/(24*7), yn(:,38), 'LineWidth', 2, 'LineStyle', linest(inh)); xlabel('Time (week)'); ylabel('Fenestration Diameter (nm)')
    z_params = params;
    end
    
    lgd = legend(['default', params{5}(inhib_prod)], 'Location', 'NorthEast');
    
    z_params = params;
    for inh = [1:length(inhib_knock_d)]
        z_params{3}(inhib_knock_d(inh)) = 0.5;
        [tn, yn] = ode15s(@coupledODE_IVV_step,tnew,y0,opts,z_params,p_params, state, GC_conc', intv); 
        yn = abs(yn);
        
        
        figure(Tstop); hold on; ax = gca; ax.FontSize = 14;
        subplot(2,2,3); hold on; plot(tn/(24*7), yn(:,38), 'LineWidth', 2, 'LineStyle', linest(inh)); xlabel('Time (week)'); ylabel('Fenestration Diameter (nm)')
        z_params = params;
    end
    
    lgd = legend(['default', params{4}(inhib_knock_d)], 'Location', 'NorthEast'); 

end
%%  4: Glucose-intervention
if task == 4
    opts=[];
    % no intervention step
    intv = 'none';
    t1 = [336:3360];
    %[T, Y] = ode15s(@coupledODE_IVV_step,t1,y0,opts,params,p_params, state, GC_conc', intv);
    %Y = real(Y);

    opts=[];
    Nn = 25;
    glu_sampled = zeros(11,Nn);

    for Nstep = 1:Nn
        for i = 1:length(GC_time)
            glu_sampled(i,Nstep) = unifrnd(GC_LB(:,i), GC_UB(:,i)); %

        end

        [T, Y] = ode15s(@coupledODE_IVV_step,t1,y0,opts,params,p_params, state, glu_sampled(:,Nstep), intv);
        YstepP(Nstep,:,:) = real(Y);
        fprintf('run %i finished\n', Nstep)

    end

    Ymean(1,:,:) = mean(YstepP([1:Nn],:,:));
    Ymean = squeeze(Ymean(1,:,:));   
       

    for Nstep=1:Nn
        %Gp = step_function(glu_sampled(:,Nstep));
        for Tt = [336:1:tspan(end)]
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
        for Tt = [336:1:tspan(end)]
            if (Tt >= time_lee(1) && Tt <= time_lee(3))
                GLU_p4(Tt,1) = 0.051*(Tt)  - 9.38;
            else
                GLU_p4(Tt,1) = 0.051*(336)  - 9.38;
            end
        end
    
    

    y0 = Y4(end,:);
    intv = '10h';
    % intervention at 10 hours 
    [T10, Y10] = ode15s(@coupledODE_IVV_step,t1,y0,opts,params,p_params, state, GC_conc', intv);
    Y10 = real(Y10);

    opts=[];
    Nn = 25;
    glu_sampled = zeros(11,Nn);

    for Nstep = 1:Nn
        for i = 1:length(GC_time)
            glu_sampled(i,Nstep) = unifrnd(GC_LB(:,i), GC_UB(:,i)); %

        end

        [T10, Y10] = ode15s(@coupledODE_IVV_step,t1,y0,opts,params,p_params, state, glu_sampled(:,Nstep), intv);
        YstepP10(Nstep,:,:) = real(Y10);
        fprintf('run %i finished\n', Nstep)

    end

    Ymean10(1,:,:) = mean(YstepP10([1:Nn],:,:));
    Ymean10 = squeeze(Ymean10(1,:,:));   
    
    for Nstep=1:Nn
        %Gp = step_function(glu_sampled(:,Nstep));
        for Tt = [336:1:tspan(end)]
            if (Tt >= time_lee(1) && Tt <= time_lee(5))
                GLU_p10(Tt,Nstep) =  0.051*(Tt)  - 9.38;
            elseif (Tt>time_lee(5) && Tt<=time_lee(9))
           
                GLU_p10(Tt,Nstep) = step_function(Tt, glu_sampled(:,Nstep));
            else
                GLU_p10(Tt,Nstep) = 0.051*(336)  - 9.38;
            end
        end
    end

    
 %%   

    figure(4); box;  
    subplot(1,2,1)
    plot(T/(24*7), GLU_p([336:3360],:), 'LineWidth', 1.5); hold on;
    hold on;
    
    xlabel('Time (weeks)'); xlim([0,21]);
    ylabel('Glucose (mmol/l)'); ylim([0,55]);
    ax = gca; ax.FontSize = 20;

    subplot(1,2,2)
    plot(T/(24*7), mean(GLU_p([336:3360],:)'), 'LineWidth', 3, 'Color', 'k'); 
    hold on; scatter(time_finch/(7*24), glu_finch, 'MarkerFaceColor',[1 0 0],'Marker','^','Color',[1 0 0])
    hold on; errorbar(time_finch/(7*24), glu_finch, abs(glu_finch - glu_LB), abs(glu_finch - glu_UB), '^', 'Color',[1  0  0]);
    hold on; scatter(time_lee/(7*24), glucose_lee, 'MarkerFaceColor',[0  0  1],'Marker','o','Color',[0  0  1])
    hold on; errorbar(time_lee/(7*24), glucose_lee, abs(glucose_lee - LB_lee), abs(glucose_lee - UB_lee), 'o', 'Color',[0  0  1]);
    
%     hold on
%     jbfill([336:3360]/24/7, min(GLU_p([336:3360],:)'), max(GLU_p([336:3360],:)'), [.7 .7 .7], [.7 .7 .7], 1, 0.2);
    hold on
    plot(T4/(24*7), GLU_p4([336:3360],1), 'LineWidth', 3, 'Color', 'k', 'LineStyle', ':'); 
    hold on;
    plot(T10/(24*7), mean(GLU_p10([336:3360],:)'), 'LineWidth', 3, 'Color', [0.7 0.7 0.7], 'LineStyle', '--'); 

    hold on;
    xlabel('Time (weeks)'); xlim([0,21]);
    ylabel('Glucose (mmol/l)'); ylim([0,55]);
    ax = gca; ax.FontSize = 20;
    hold on
    legend('Model (no intervention)', 'Finch et al. (2022)', '', 'Lee et al. (2018)', '',  'intervention at 4 weeks', 'intervention at 10 weeks');
    ylim([0,52]);
    x0=10;
    y0=10;
    width=800;
    height=800;
    set(gcf,'position',[x0,y0,width,height])
    ax = gca; ax.FontSize = 20;

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
    subplot(3,2,1);
    bar((Ymean(end,[2:29,31:36]) - Ymean(1,[2:29,31:36])), 'k'); title('No intervention'); 
    xticks([1:34]); xticklabels(params{4}([2:29,31:36])); ylabel('Change relative to baseline')
    
    figure(5)
    subplot(3,2,2)
    bar((Ymean(end,[37:38]) - Ymean(1,[37:38])), 'k'); title('No intervention')
    xticks([1:2]); xticklabels(params{4}([37:38])); ylabel('Change relative to baseline')

       hold on
    
    figure(5);
    subplot(3,2,3);
    bar((Y4(end,[2:29,31:36]) - Y4(1,[2:29,31:36])), 'k'); title('Glucose intervention at 4 weeks'); 
    xticks([1:34]); xticklabels(params{4}([2:29,31:36])); ylabel('Change relative to baseline')
    
    figure(5)
    subplot(3,2,4)
    bar((Y4(end,[37:38]) - Y4(1,[37:38])), 'k'); title('Glucose intervention at 4 weeks')
    xticks([1:2]); xticklabels(params{4}([37:38])); ylabel('Change relative to baseline')

       hold on
    
    figure(5);
    subplot(3,2,5);
    bar((Ymean10(end,[2:29,31:36]) - Ymean10(1,[2:29,31:36])), 'k'); title('Glucose intervention at 10 weeks'); 
    xticks([1:34]); xticklabels(params{4}([2:29,31:36])); ylabel('Change relative to baseline')
    
    figure(5)
    subplot(3,2,6)
    bar((Ymean10(end,[37:38]) - Ymean10(1,[37:38])), 'k'); title('Glucose intervention at 10 weeks')
    xticks([1:2]); xticklabels(params{4}([37:38])); ylabel('Change relative to baseline')
end