% Author: Krutika Patidar
% Dated: August 30th 2021
% Description: The function evaluates minimum root mean squared error

function [min_error, weighted_err]= coupledODE_error(p, params, y0, tspan, tau_index, n_index, k_index)



% Uncomment if parameters (p) are being optimized

% size_tau = size(tau_index,2);
% size_n = size(n_index,2);
% size_k = size(k_index,2);

%for m = 1:size_tau
%        params{2}(tau_index(m)) = p(m)*(max(params{2})-min(params{2})) + min(params{2});       
%end
%for n = 1:size_n
%        params{1}(2,n_index(n)) = p(n + size_tau)*(max(params{1}(2,:)) -  min(params{1}(2,:))) + min(params{1}(2,:));
%end
%for o = 1:size_k % EC50
%        params{1}(3,k_index(o)) = p(o+size_tau + size_n);
%end


options = [];


[t, y] = ode23(@coupledODE_opt_extended,tspan,y0,options,params);


Yout = y;

load invitro_data.mat

% Model training and validation
 GLU = params{1}(1,1);
 LPS = params{1}(1,2);
 
 error = zeros(1,8);
 w = zeros(1,8);


 if GLU==1 && LPS==0 
         
         error(1,1) = sum((Yout([1,6,24], 23) - train_data([2,6,10],1)).^2); 
         error(1,2) = sum((Yout([1,48], 25) - train_data([2,14],2)).^2);
         error(1,3) = sum((Yout([1,6,24], 24) - train_data([2,6,10],3)).^2); 
         error(1,4) = sum((Yout([1,48],6) - train_data([2,14],4)).^2) ; 
         error(1,5) = sum((Yout([1,24], 13) - train_data([2,14],5)).^2); 
         error(1,6) = sum((Yout([1,24], 22) - train_data([2,10],6)).^2); 
         error(1,7) = sum((Yout([1,3,6,12,24],20) - train_data([2,6,10,14,18],7)).^2); 
         error(1,8) = sum((Yout([1,12,24,48],12) - train_data([2,6,10,14],9)).^2); 
         

         % weight choices: 1/Y_out, 1/Yout^2, 1/train_data,
         % 1/train_data^2

          w(1,1) = 1/std(train_data([2,6,10],1));
          w(1,2) = 1/std(train_data([2,14],2));
          w(1,3) = 1/std(train_data([2,6,10],3));
          w(1,4) = 1/std(train_data([2,14],4));
          w(1,5) = 1/std(train_data([2,14],5));
          w(1,6) = 1/std(train_data([2,10],6));
          w(1,7) = 1/std(train_data([2,6,10,14,18],7));
          w(1,8) = 1/std(train_data([2,6,10,14],9));
 end

 if GLU==0 && LPS==1 
         error(1,1) = sum((Yout([1,6,24], 23) - train_data([3,7,11],1)).^2);     %IL-6 
         error(1,2) = sum((Yout([1,48], 25) - train_data([3,15],2)).^2);         %IL-1b
         error(1,3) = sum((Yout([1,6,24], 24) - train_data([3,7,11],3)).^2);     %TNF-a
         error(1,4) = sum((Yout([1,48],6) - train_data([3,15],4)).^2);           %VEGFa-mRNA
         error(1,5) = sum((Yout([1,24], 13) - train_data([3,15],5)).^2);         %ROS 
         error(1,6) = sum((Yout([1,24], 22) - train_data([3,11],6)).^2);         %eNOS
         error(1,7) = sum((Yout([1,48],20) - train_data([3,23],7)).^2);          % NO 
         error(1,8) = sum((Yout([1,48],12) - train_data([3,15],9)).^2);          %ROSec
        
         
          w(1,1) = 1/std(train_data([3,7,11],1));
          w(1,2) = 1/std(train_data([3,15],2));
          w(1,3) = 1/std(train_data([3,7,11],3));
          w(1,4) = 1/std(train_data([3,15],4));
          w(1,5) = 1/std(train_data([3,15],5));
          w(1,6) = 1/std(train_data([3,11],6));
          w(1,7) = 1/std(train_data([3,23],7));
          w(1,8) = 1/std(train_data([3,15],9));
         
 end
 if  GLU==1 && LPS==1 
         error(1,1) = sum((Yout([1,6,24], 23) - train_data([4,8,12],1)).^2); 
         error(1,2) = sum((Yout([1,48], 25) - train_data([4,16],2)).^2); 
         error(1,3) = sum((Yout([1,6,24], 24) - train_data([4,8,12],3)).^2); 
         error(1,4) = sum((Yout([1,48],6) - train_data([4,16],4)).^2);
         error(1,5) = sum((Yout([1,24], 13) - train_data([4,16],5)).^2); 
         error(1,6) = sum((Yout([1,24], 22) - train_data([4,12],6)).^2); 
         error(1,7) = sum((Yout([1,48],20) - train_data([4,24],7)).^2) ; 
         error(1,8) = sum((Yout([1,12,24,48],12) - train_data([4,8,12,16],9)).^2);  
         
          w(1,1) = 1/std(train_data([4,8,12],1));
          w(1,2) = 1/std(train_data([4,16],2));
          w(1,3) = 1/std(train_data([4,8,12],3));
          w(1,4) = 1/std(train_data([4,16],4));
          w(1,5) = 1/std(train_data([4,16],5));
          w(1,6) = 1/std(train_data([4,12],6));
          w(1,7) = 1/std(train_data([4,24],7));
          w(1,8) = 1/std(train_data([4,8,12,16],9));
 end
  
  
 min_error = (error(1,1) + error(1,2) + error(1,3) + error(1,4) + error(1,5) + error(1,6) + error(1,7)); 
 weighted_err = w(1)*(error(1,1)) + w(2)*(error(1,2)) + w(3)*(error(1,3)) + w(4)*(error(1,4)) + w(5)*(error(1,5)) + w(6)*(error(1,6)) + w(7)*(error(1,7)) + w(8)*error(1,8);
 


end