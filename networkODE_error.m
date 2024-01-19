% Author: Krutika Patidar
% Dated: August 30th 2021
% Description: The function evaluates minimum sum of squared error for the network model

function [min_error] = networkODE_error(p, params, y0, tspan, tau_index, n_index, k_index, W_index)

size_tau = size(tau_index,2);
size_n = size(n_index,2);
size_k = size(k_index,2);
size_W = size(W_index,2);


for tau_iter = 1:size_tau
        params{2}(tau_index(tau_iter)) = p(tau_iter);      
end
for w_iter = 1:size_W 
        params{1}(1, W_index(w_iter)) = p(w_iter +  size_tau);
end
for n_iter = 1:size_n
        params{1}(2,n_index(n_iter)) = p(n_iter + size_tau + size_W); 
end
for k_iter = 1:size_k 
        params{1}(3,k_index(k_iter)) = p(k_iter + size_tau + size_n + size_W);
end

options = [];


[t, y] = ode23s(@networkODE,tspan,y0,options,params);


Yout = real(y); Tout = t;


load invitro_data.mat

% Model training and validation
 GLU = params{1}(1,1);
 LPS = params{1}(1,2);
 error = zeros(1,8);
 pred = zeros(1,8);


 if GLU==1 && LPS==0 % Condition (1,0)
            
        
         %error(1,1) = sum((Yout([1,6,24], 23) - train_data([2,6,10],1)).^2); 
         %error(1,2) = sum((Yout([1,48], 25) - train_data([2,14],2)).^2);
         %error(1,3) = sum((Yout([1,6,24], 24) - train_data([2,6,10],3)).^2); 
         %error(1,4) = sum((Yout([1,48],6) - train_data([2,14],4)).^2) ; 
         %error(1,5) = sum((Yout([1,24], 13) - train_data([2,14],5)).^2); 
         %error(1,6) = sum((Yout([1,24], 22) - train_data([2,10],6)).^2); 
         %error(1,7) = sum((Yout([1,3,6,12,24],20) - train_data([2,6,10,14,18],7)).^2); 
         %error(1,8) = sum((Yout([1,12,24,48],12) - train_data([2,6,10,14],9)).^2); 
         
         pred_in = interp1(Tout, Yout(:,23), time_data([2,6,10],1)); error(1,1) = sum((pred_in - train_data([2,6,10],1)).^2); pred_in = [];
         pred_in = interp1(Tout, Yout(:,25), time_data([2,14],2)); error(1,2) = sum((Yout([1,48], 25) - train_data([2,14],2)).^2); pred_in = [];
         pred_in = interp1(Tout, Yout(:,24), time_data([2,6,10],3)); error(1,3) = sum((Yout([1,6,24], 24) - train_data([2,6,10],3)).^2); pred_in = [];
         pred_in = interp1(Tout, Yout(:,6), time_data([2,14],4)); error(1,4) = sum((Yout([1,48],6) - train_data([2,14],4)).^2) ; pred_in = [];
         pred_in = interp1(Tout, Yout(:,13), time_data([2,14],5)); error(1,5) = sum((Yout([1,24], 13) - train_data([2,14],5)).^2); pred_in = [];
         pred_in = interp1(Tout, Yout(:,22), time_data([2,10],6)); error(1,6) = sum((Yout([1,24], 22) - train_data([2,10],6)).^2); pred_in = [];
         pred_in = interp1(Tout, Yout(:,20), time_data([2,6,10,14,18],7)); error(1,7) = sum((Yout([1,3,6,12,24],20) - train_data([2,6,10,14,18],7)).^2);pred_in = []; 
         pred_in = interp1(Tout, Yout(:,12), time_data([2,6,10,14],9)); error(1,8) = sum((Yout([1,12,24,48],12) - train_data([2,6,10,14],9)).^2); pred_in = [];
         
          
 end

 if GLU==0 && LPS==1 % Condition (0,1)
         
         pred_in = interp1(Tout, Yout(:,23), time_data([3,7,11],1)); error(1,1) = sum((Yout([1,6,24], 23) - train_data([3,7,11],1)).^2); pred_in = [];    %IL-6 
         pred_in = interp1(Tout, Yout(:,25), time_data([3,15],2)); error(1,2) = sum((Yout([1,48], 25) - train_data([3,15],2)).^2);   pred_in = [];      %IL-1b
         pred_in = interp1(Tout, Yout(:,24), time_data([3,7,11],3)); error(1,3) = sum((Yout([1,6,24], 24) - train_data([3,7,11],3)).^2); pred_in = [];    %TNF-a
         pred_in = interp1(Tout, Yout(:,6), time_data([3,15],4)); error(1,4) = sum((Yout([1,48],6) - train_data([3,15],4)).^2);   pred_in = [];         %VEGFa-mRNA
         pred_in = interp1(Tout, Yout(:,13), time_data([3,15],5)); error(1,5) = sum((Yout([1,24], 13) - train_data([3,15],5)).^2);  pred_in = [];       %ROS 
         pred_in = interp1(Tout, Yout(:,22), time_data([3,11],6)); error(1,6) = sum((Yout([1,24], 22) - train_data([3,11],6)).^2); pred_in = [];        %eNOS
         pred_in = interp1(Tout, Yout(:,20), time_data([3,23],7)); error(1,7) = sum((Yout([1,48],20) - train_data([3,23],7)).^2);  pred_in = [];        % NO 
         pred_in = interp1(Tout, Yout(:,12), time_data([3,15],9)); error(1,8) = sum((Yout([1,48],12) - train_data([3,15],9)).^2);  pred_in = [];        %ROSec
          
             
         
 end

 if  GLU==1 && LPS==1 %Condition (1,1)
         pred_in = interp1(Tout, Yout(:,23), time_data([4,8,12],1)); error(1,1) = sum((Yout([1,6,24], 23) - train_data([4,8,12],1)).^2); pred_in = [];
         pred_in = interp1(Tout, Yout(:,25), time_data([4,16],2)); error(1,2) = sum((Yout([1,48], 25) - train_data([4,16],2)).^2); pred_in = [];
         pred_in = interp1(Tout, Yout(:,24), time_data([4,8,12],3)); error(1,3) = sum((Yout([1,6,24], 24) - train_data([4,8,12],3)).^2); pred_in = [];
         pred_in = interp1(Tout, Yout(:,6), time_data([4,16],4)); error(1,4) = sum((Yout([1,48],6) - train_data([4,16],4)).^2); pred_in = [];
         pred_in = interp1(Tout, Yout(:,13), time_data([4,16],5)); error(1,5) = sum((Yout([1,24], 13) - train_data([4,16],5)).^2); pred_in = [];
         pred_in = interp1(Tout, Yout(:,22), time_data([4,12],6)); error(1,6) = sum((Yout([1,24], 22) - train_data([4,12],6)).^2); pred_in = [];
         pred_in = interp1(Tout, Yout(:,20), time_data([4,24],7)); error(1,7) = sum((Yout([1,48],20) - train_data([4,24],7)).^2) ; pred_in = [];
         pred_in = interp1(Tout, Yout(:,12), time_data([4,8,12,16],9)); error(1,8) = sum((Yout([1,12,24,48],12) - train_data([4,8,12,16],9)).^2);  pred_in = [];
       
         
               
 end
 
  
 min_error = sum(error); 



end