function [p_params] = coupledODE_physParams(state)

if state == "norm_mice"



elseif state == "diab_mice"


%     yss = 7;
%     ns = 3.99;
%     kform = 1.013;
%     kd = 65.9; 
%     ke = 2.04;
%     kloss = 4.6;
%     yss2 = 4.07;

    yss = 7; % unitless
    ns = 3.99; % unitless
    kform = 1.013; % 1/hr
    kd = 65.9; % nm/hr
    ke = 2.04; % % 1/hr
    kloss = 4.6; % 1/hr
    yss2 = 4.02; % unitless

end

p_params = [yss, ns, kform, kd, ke, kloss, yss2]; 

end