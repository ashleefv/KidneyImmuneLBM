function [params,y0] = networkODE_opt_loadParams(choice) 

% species parameters 
speciesNames = {'GLU','LPS','AGE','VEGFR1','VEGFR2','VEGF-A_m_R_N_A','RAGE_e_c','RAGE','TLR4','NADPH','NADPH_e_c','ROS_e_c','ROS','PI3K','AKT','PI3K_e_c','AKT_e_c','NF\kappaB_e_c','NF\kappaB','NO','ONOO','eNOS','IL-6','TNF-\alpha','IL-1\beta','PLC-\gamma','VEGF-A','pJunction','Ca','Gap Width',}; 
tau = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ]; 
ymax = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ];
 
% reaction parameters optimized from multistart_param_opt.m (means of acceptable paramters from GLU only condition)
w = [0	0	1	1	1	0.944	1	0.950	0.943	1	0.950	0.999	1	0.950	0.950	0.950	0.949	0.950	1	1	1	0.938	1	1	0.946	1	1	1	1	1	0.986	0.962	1	1	1	1	1	1	1	0.01];
n = [1.4	1.4	2.71	1.45	2.71	2.70	2.70	2.70	2.64	2.70	2.70	3.96	2.70	2.70	2.70	2.70	2.77	2.70	2.71	2.72	1.56	1.60	2.73	2.72	3.91	2.93	3.11	1.4	1.4	2.70	2.70	2.70	1.46	3.66	1.4	2.68	1.4	1.4	1.4	1.4];
EC50 = [0.5	0.5	0.419	0.5	0.474	0.470	0.5	0.422	0.545	0.5	0.419	0.840	0.5	0.420	0.419	0.422	0.595	0.421	0.5	0.5	0.839	0.839	0.5	0.5	0.838	0.688	0.576	0.5	0.5	0.010	0.471	0.391	0.157	0.833	0.5	0.5	0.5	0.5	0.5	0.5];

y0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ]; 

% OPTIMIZED SPECIES PARAMS (LISTED BELOW) optimized from multistart_param_opt.m (means of acceptable paramters from each condition)
if choice=="GLU"
    disp("optimal species params: GLU");
    % GLUCOSE ONLY
    y0(27) = 1; y0(6) = 1; y0(20) = 1; y0(22) = 1;
    tau(:) = [1	5.06	1	1	1	3.39	1	1	5.04	1	1	9.97	5.89	5.04	1	6.60	5.75	1	5.07	4.36	1	1.91	4.66	9.90	4.90	1	3.14	1	1	1];

end
if choice=="LPS"
    % LPS ONLY
    disp("optimal species params: LPS");
    tau(:) = [1	0.800	1	1	1	0.944	1	1	0.800	1	1	1.07	0.800	0.800	1	0.977	0.984	1	0.800	1.00	1	0.981	0.800	0.800	1.00	1	0.971	1	1	1];
end

if choice=="both"
    % GLUCOSE + LPS ONLY
    disp("optimal species params: both");
    tau(:) = [1	1.00	1	1	1	1.48	1	1	1.00	1	1	5.15	1.00	1.00	1	1.00	1.00	1	1.00	6.80	1	1.00	1.00	1.00	4.95	1	1.10	1	1	1];
end


 
rpar = [w;n;EC50];

params = {rpar,tau,ymax,speciesNames};
