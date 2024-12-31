function [params,y0] = coupledODE_IVIVC_params() 
% species parameters 
speciesNames = {'GLU','Actin_s','AGE','VEGFR1','VEGFR2','VEGF-A_{mRNA}','RAGE_{ec}','RAGE','IL-1R','NADPH','NADPH_{ec}','ROS_{ec}','ROS','PI3K','AKT','PI3K_{ec}','AKT_{ec}','NF\kappaB_{ec}','NF\kappaB','NO','ONOO','eNOS','IL-6','TNF-\alpha','IL-1\beta','PLC-\gamma','VEGF-A','pJunction','Ca','GapWidth','Actin_r','RhoRock','MLCK','pMLC','MLCP','MLC','Number','Diameter'}; 

tau = [1	1	1	0.35	0.35	88	0.35	0.35	0.35	1	1	1	1	1	1	1	1	0.055	0.055	1	1	1	90	90	90	1	1.13  1	 1	1  1  1  1  1  1  1  1  1]; 
ymax = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

reac_names = {'\Rightarrow GLU','IL-1\beta \Rightarrow IL-1R','GLU \Rightarrow AGE', 'AGE \Rightarrow RAGE','RAGE \Rightarrow NADPH', 'IL-1R & ROS \Rightarrow NF\kappaB', 'IL-1R \Rightarrow PI3K','NADPH \Rightarrow ROS', 'PI3K \Rightarrow AKT', 'PI3K \Rightarrow ROS', 'NF\kappaB_{ec} \Rightarrow TNF-\alpha','AKT \Rightarrow NF\kappaB','NF\kappaB \Rightarrow IL-6','NF\kappaB \Rightarrow TNF-\alpha','NF\kappaB \Rightarrow VEGF-A_{mRNA}', 'VEGF-A_{mRNA} \Rightarrow VEGF-A','NF\kappaB \Rightarrow IL-1\beta','VEGF-A \Rightarrow VEGFR1','VEGF-A \Rightarrow VEGFR2','AGE \Rightarrow RAGE_{ec}','RAGE_{ec} \Rightarrow NADPH_{ec}','VEGFR2 \Rightarrow PI3K_{ec}','VEGFR1 \Rightarrow PI3K_{ec}', 'NADPH_{ec} \Rightarrow ROS_{ec}', 'PI3K_{ec} \Rightarrow AKT_{ec}', 'AKT_{ec} \Rightarrow eNOS' , 'VEGFR1 \Rightarrow PLC-\gamma' ,'PLC-\gamma \Rightarrow NF\kappaB_{ec}' , 'ROS_{ec} \Rightarrow NF\kappaB_{ec}','NF\kappaB_{ec} \Rightarrow IL-6', 'NF\kappaB_{ec} \Rightarrow IL-1\beta', 'eNOS  \Rightarrow NO', 'eNOS \Rightarrow ROS_{ec}' ,'ROS_{ec} & NO \Rightarrow ONOO','!NO \Rightarrow Ca','PLC-\gamma =\Rightarrow Ca', 'Ca \Rightarrow pJunction','pJunction \Rightarrow Gap Width', 'Ca \Rightarrow NO', '!NO & ROS \Rightarrow MLCK', 'Ca \Rightarrow MLCK', 'RhoRock \Rightarrow pMLC', 'pMLC & MLCP \Rightarrow MLC', 'VEGFR2 \Rightarrow RhoRock','!RhoRock \Rightarrow MLCP', 'MLC & MLCK \Rightarrow pMLC', 'pMLC \Rightarrow Actin_s', 'MLC \Rightarrow Actin_r'}; % 'pMLC => Fen Diameter', 'ActinR & !ActinS \Rightarrow Fenestration Count'};
 
% reaction parameters optimized from multistart_param_opt.m (means of acceptable paramters from GLU only condition)
w = [0	1	1	1	0.944	1	1	0.943	1	0.950	0.999	1	0.950	0.950	0.950	0.949	0.950	1	1	1	0.938	1	1	0.946	1	1	1	1	1	0.986	0.962	1  1  1  1  1  1  1  1 1 1 1 1 1 1 1 1 1];
n = [1.4  1.4	1.45	2.71	2.70	1.4	1.4	2.64	2.70	2.70	3.96	2.70	2.70	2.70	2.70	2.77	2.70	2.71	2.72	1.56	1.60	2.73	2.72	3.91	2.93	3.11	1.4	1.4	2.70	2.70	2.70	1.46	3.66	1.4	2.68	1.4	1.4	1.4	1.4   1.4  1.4  1.4 1.4	1.4	1.4	1.4   1.4  1.4];
EC50 = [0.5 0.5	0.5	0.474	0.470	0.5	0.5	0.545	0.5	0.419	0.840	0.5	0.420	0.419	0.422	0.595	0.421	0.5	0.5	0.839	0.839	0.5	0.5	0.838	0.688	0.576	0.5	0.5	0.010	0.471	0.391	0.157	0.833	0.5	0.5	0.5	0.5	0.5	0.5  0.5  0.5  0.5  0.5	0.5	0.5	0.5	0.5	0.5];


y0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; 

load data/CGM_db.mat

 
y0(37) = 6.3;
y0(38) = 47.91;
tau([34]) = [400.5];
 

rpar = [w;n;EC50];



params = {rpar,tau,ymax,speciesNames, reac_names};
end
