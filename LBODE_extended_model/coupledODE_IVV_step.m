function [dydt] = coupledODE_IVV_step(t,y,params,p_params,state, glu_sampled, intv) 

% state = "diab_mice"; % note GLU params work better for T2D_human state

global number_ctrl time_ctrl density diameter GC_conc GC_time GC_LB GC_UB time_lee glu_UB glu_LB time_finch glu_finch LB_lee UB_lee glucose_lee

GLU = 1; 
ActinS = 2; %
AGE = 3; 
VEGFR1 = 4; 
VEGFR2 = 5; 
VEGFamRNA = 6; 
RAGEec = 7; 
RAGE = 8; 
IL1R = 9; 
NADPH = 10; 
NADPHec = 11; 
ROSec = 12; 
ROS = 13; 
PI3K = 14; 
AKT = 15; 
PI3Kec = 16; 
AKTec = 17; 
NFKBec = 18; 
NFKB = 19; 
NO = 20; 
ONOO = 21; 
eNOS = 22; 
IL6 = 23; 
TNFa = 24; 
IL1b = 25; 
PLC = 26; 
VEGFa = 27; 
pJunc = 28; 
Calcium = 29; 
%GapWidth = 30; 
ActinR = 30;
RhoRock = 31;
MLCK = 32;
pMLC = 33;
MLCP = 34;
MLC = 35;
FenCount = 36;
FenDiameter = 37; 

dydt = zeros(37,1); 

rpar = params{1};
tau = params{2};
ymax = params{3};




if state == 'diab_mice'
    Gp0 =   0.051*(336)  - 9.38;
    %Gp = step_function(glu_sampled);
if intv == "none"
    if (t>=time_lee(1) && t<=time_lee(5))
       Gp1 =   0.051*(t) - 9.38;
       W_GLU = (Gp1 - Gp0) / (max(glu_UB) - Gp0);  % (1) use of max of finch et al. data
    else
       Gp = step_function(t, glu_sampled);
       %W_GLU = (double(Gp([t])) - Gp0) / (max(GC_conc) - Gp0);  
       W_GLU = (Gp - Gp0) / (max(glu_UB) - Gp0);   
    end
elseif intv == "4h"

    if (t>=time_lee(1) && t<=time_lee(3))
       Gp1 =   0.051*(t) - 9.38;
       W_GLU = (Gp1 - Gp0) / (max(glu_UB) - Gp0);  % (1) use of max of finch et al. data
    else
        Gp1 = Gp0;
        W_GLU = (Gp1 - Gp0) / (max(glu_UB) - Gp0);  % (1) use of max of finch et al. data
    end

elseif intv == "10h"
   
    if (t>=time_lee(1) && t<=time_lee(5))
       Gp1 =   0.051*(t) - 9.38;
       W_GLU = (Gp1 - Gp0) / (max(glu_UB) - Gp0);  % (1) use of max of finch et al. data
    elseif (t>time_lee(5) && t<=time_lee(9))

        Gp1 = step_function(t, glu_sampled);
        W_GLU = (Gp1 - Gp0) / (max(glu_UB) - Gp0);  % (1) use of max of finch et al. data
    else 
        Gp1 = Gp0;
        W_GLU = (Gp1 - Gp0) / (max(glu_UB) - Gp0); 
    end
   
    
end
end

if state == 'norm_mice'
     Gp0 =   0.051*(336)  - 9.38;
     Gp1 = Gp0;
     W_GLU = (Gp1 - Gp0) / (max(glu_UB) - Gp0);
     
end

dydt(GLU) = (ymax(GLU)*W_GLU - y(GLU))/tau(GLU);   % LBM glucose 

% time-dependent glucose
%
dydt(AGE) = (act(y(GLU),rpar(:,3))*ymax(AGE) - y(AGE))/tau(AGE); 
dydt(VEGFR1) = (act(y(VEGFa),rpar(:,18))*ymax(VEGFR1) - y(VEGFR1))/tau(VEGFR1); 
dydt(VEGFR2) = (act(y(VEGFa),rpar(:,19))*ymax(VEGFR2) - y(VEGFR2))/tau(VEGFR2); 
dydt(VEGFamRNA) = (act(y(NFKB),rpar(:,17))*ymax(VEGFamRNA) - y(VEGFamRNA))/tau(VEGFamRNA); 
dydt(RAGEec) = (act(y(AGE),rpar(:,20))*ymax(RAGEec) - y(RAGEec))/tau(RAGEec); 
dydt(RAGE) = (act(y(AGE),rpar(:,4))*ymax(RAGE) - y(RAGE))/tau(RAGE); 

dydt(IL1R) = (act(y(IL1b),rpar(:,2))*ymax(IL1R) - y(IL1R))/tau(IL1R);

dydt(NADPH) = (act(y(RAGE),rpar(:,5))*ymax(NADPH) - y(NADPH))/tau(NADPH); 
dydt(NADPHec) = (act(y(RAGEec),rpar(:,21))*ymax(NADPHec) - y(NADPHec))/tau(NADPHec); 

dydt(ROSec) = (OR(act(y(NADPHec),rpar(:,24)),act(y(eNOS),rpar(:,34)))*ymax(ROSec) - y(ROSec))/tau(ROSec); 
dydt(ROS) = (OR(act(y(NADPH),rpar(:,8)),act(y(PI3K),rpar(:,10)))*ymax(ROS) - y(ROS))/tau(ROS); 
dydt(PI3K) = (act(y(IL1R),rpar(:,7))*ymax(PI3K) - y(PI3K))/tau(PI3K); 

%dydt(ROS) = (act(y(NADPH),rpar(:,9))*ymax(ROS) - y(ROS))/tau(ROS); % here
%dydt(PI3K) = (OR(act(y(ROS),rpar(:,11)), act(y(TLR4),rpar(:,8)))*ymax(PI3K) - y(PI3K))/tau(PI3K); % here

dydt(AKT) = (act(y(PI3K),rpar(:,9))*ymax(AKT) - y(AKT))/tau(AKT); 
dydt(PI3Kec) = (OR(act(y(VEGFR2),rpar(:,22)),act(y(VEGFR1),rpar(:,23)))*ymax(PI3Kec) - y(PI3Kec))/tau(PI3Kec); 
dydt(AKTec) = (act(y(PI3Kec),rpar(:,25))*ymax(AKTec) - y(AKTec))/tau(AKTec); 
dydt(NFKBec) = (OR(act(y(PLC),rpar(:,28)),act(y(ROSec),rpar(:,29)))*ymax(NFKBec) - y(NFKBec))/tau(NFKBec); 

dydt(NFKB) = (OR(AND(rpar(:,6),act(y(IL1R),rpar(:,6)),act(y(ROS),rpar(:,6))),act(y(AKT),rpar(:,12)))*ymax(NFKB) - y(NFKB))/tau(NFKB); 
%dydt(NFKB) = (OR(act(y(ROS),rpar(:,7)),act(y(AKT),rpar(:,13)))*ymax(NFKB) - y(NFKB))/tau(NFKB); 


dydt(NO) = (OR(act(y(eNOS),rpar(:,31)),act(y(Calcium),rpar(:,38)))*ymax(NO) - y(NO))/tau(NO); %%%(AND(rpar(:,33),act(y(eNOS),rpar(:,33)),act(y(Calcium),rpar(:,33)))*ymax(NO) - y(NO))/tau(NO); % 
dydt(ONOO) = (AND(rpar(:,34),act(y(ROSec),rpar(:,34)),act(y(NO),rpar(:,34)))*ymax(ONOO) - y(ONOO))/tau(ONOO); 
% dydt(NOu) = (y(NO) - AND(rpar(:,35),act(y(ROSec),rpar(:,35)),act(y(NO),rpar(:,35)))*ymax(NOu))/(taus(NOu)) - (y(NOu))/tau(NOu) - (y(NO)*y(ROSec))/tau(NOu) ; 

dydt(eNOS) = (act(y(AKTec),rpar(:,26))*ymax(eNOS) - y(eNOS))/tau(eNOS); 

dydt(IL6) = (OR(act(y(NFKBec),rpar(:,30)), act(y(NFKB),rpar(:,13)))*ymax(IL6) - y(IL6))/tau(IL6); % or t/b
dydt(TNFa) = (OR(act(y(NFKBec),rpar(:,11)),act(y(NFKB),rpar(:,14)))*ymax(TNFa) - y(TNFa))/tau(TNFa); 

dydt(IL1b) = (OR(act(y(NFKBec),rpar(:,31)), act(y(NFKB),rpar(:,17)))*ymax(IL1b) - y(IL1b))/tau(IL1b); 
dydt(PLC) = (act(y(VEGFR1),rpar(:,27))*ymax(PLC) - y(PLC))/tau(PLC); 
dydt(VEGFa) = (act(y(VEGFamRNA),rpar(:,16))*ymax(VEGFa) - y(VEGFa))/tau(VEGFa);
dydt(Calcium) = (OR(inhib(y(NO),rpar(:,35)),act(y(PLC),rpar(:,36)))*ymax(Calcium) - y(Calcium))/tau(Calcium); % calcium oscillations follow glucose oscillations
  


dydt(pJunc) = (act(y(Calcium),rpar(:,37))*ymax(pJunc) - y(pJunc))/tau(pJunc);
dydt(MLCK) = (OR(AND(rpar(:,39),act(y(ROS),rpar(:,39)),inhib(y(NO),rpar(:,39))), act(y(Calcium),rpar(:,40)))*ymax(MLCK) - y(MLCK))/tau(MLCK); 
dydt(pMLC) = (OR(act(y(RhoRock),rpar(:,41)),AND(rpar(:,45),act(y(MLCK),rpar(:,45)),act(y(MLC),rpar(:,45))))*ymax(pMLC) - y(pMLC))/tau(pMLC); 

dydt(MLCP) = (inhib(y(RhoRock),rpar(:,44))*ymax(MLCP) - y(MLCP))/tau(MLCP); 
% dydt(MLCP) = (act(y(RhoRock),rpar(:,44))*ymax(MLCP) - y(MLCP))/tau(MLCP); % healthy state simulation for Fenestration Formation

dydt(MLC) = (AND(rpar(:,42),act(y(pMLC),rpar(:,42)),act(y(MLCP),rpar(:,42)))*ymax(MLC) - y(MLC))/tau(MLC); 
dydt(RhoRock) = (act(y(VEGFR2),rpar(:,43))*ymax(RhoRock) - y(RhoRock))/tau(RhoRock); 



dydt(ActinS) = (act(y(pMLC),rpar(:,46))*ymax(ActinS) - y(ActinS))/tau(ActinS); 
dydt(ActinR) = (act(y(MLC),rpar(:,47))*ymax(ActinR) - y(ActinR))/tau(ActinR); 


dydt(FenCount) = p_params(3)*y(ActinR)*(abs(p_params(1) - y(FenCount)))^p_params(2) - p_params(6)*(y(ActinS))*(abs(p_params(7) - y(FenCount)))^p_params(2);

y0 = 47.91;
dydt(FenDiameter) = p_params(4)*(y(pMLC) - 0)^p_params(2) - p_params(5)*(y(FenDiameter) - y0);



end

% utility functions 
function fact = act(x,rpar) 
% hill activation function with parameters w (weight), n (Hill coeff), EC50 
    w = rpar(1); 
    n = rpar(2); 
    EC50 = rpar(3); 
    beta = (EC50.^n - 1)./(2*EC50.^n - 1); 
    K = (beta - 1).^(1./n); 
    fact = w.*(beta.*x.^n)./(K.^n + x.^n); 
    if fact> w                 % cap fact(x)<= 1 
        fact = w; 
    end
    
end

function finhib = inhib(x,rpar) 
% inverse hill function with parameters w (weight), n (Hill coeff), EC50 
    finhib = rpar(1) - act(x,rpar);
end

function z = OR(x,y) 
% OR logic gate 
    z = x + y - x*y;
end

function z = AND(rpar,varargin) 
% AND logic gate, multiplying all of the reactants together 
    w = rpar(1); 
    if w == 0 
        z = 0; 
    else 
        v = cell2mat(varargin); 
        z = prod(v)/w^(nargin-2);  
    end 
end
