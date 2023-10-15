function Y = UQLab_network_run(X)

firstTime = 0;
lastTime = 48;    % Duration time of simulation.

nSteps = numel(48); % Number of timesteps = 1,  Number of timesteps = 2 when [0:48:48]
%% Initialize
nReal = size(X,1); 
disp(nReal)

% Parameters reference
    tau = X(:,[1:30]);
    w = X(:,[31:70]);
    k = X(:,[71:110]);
    n = X(:,[111:150]);

% known parameters
ymax = 1;


% Initial conditions
y0 = zeros(30,1);

%% Solve equation

% solver options (smaller tolerance)
odeOpts = odeset('RelTol',1e-4,'AbsTol',1e-7);
%odeOpts = odeset('RelTol',2e-3,'AbsTol',2e-6);

% for loop to solve equations with multiple initial values and parameters
Y = zeros(nReal,8*nSteps); % 1 indicates no. of component (variable)
for ii = 1:nReal
   % setup diff equations 
    diffEq = @(t,y) [  (w(ii,1)*ymax - y(1))/tau(ii,1); 
                     (w(ii,2)*ymax - y(2))/tau(ii,2);
                     ((w(ii,4)*(k(ii,4)^n(ii,4)+1)*y(1)^n(ii,4)*ymax)/(k(ii,4)^n(ii,4) + y(1)^n(ii,4)) - y(3))/tau(ii,3);
                     ((w(ii,19)*(k(ii,19)^n(ii,19)+1)*y(27)^n(ii,19)*ymax)/(k(ii,19)^n(ii,19) + y(27)^n(ii,19))- y(4))/tau(ii,4); 
                     ((w(ii,20)*(k(ii,20)^n(ii,20)+1)*y(27)^n(ii,20)*ymax)/(k(ii,20)^n(ii,20) + y(27)^n(ii,20)) - y(5))/tau(ii,5);
                     ((w(ii,16)*(k(ii,16)^n(ii,16)+1)*y(19)^n(ii,16)*ymax)/(k(ii,16)^n(ii,16) + y(19)^n(ii,16)) - y(6))/tau(ii,6); 
                     ((w(ii,21)*(k(ii,21)^n(ii,21)+1)*y(3)^n(ii,21)*ymax)/(k(ii,21)^n(ii,21) + y(3)^n(ii,21)) - y(7))/tau(ii,7); 
                     ((w(ii,5)*(k(ii,5)^n(ii,5)+1)*y(3)^n(ii,5)*ymax)/(k(ii,5)^n(ii,5) + y(3)^n(ii,5)) - y(8))/tau(ii,8); 
                     ((w(ii,3)*(k(ii,3)^n(ii,3)+1)*y(2)^n(ii,3)*ymax)/(k(ii,3)^n(ii,3) + y(2)^n(ii,3)) - y(9))/tau(ii,9); 
                     ((w(ii,6)*(k(ii,6)^n(ii,6)+1)*y(8)^n(ii,6)*ymax)/(k(ii,6)^n(ii,6) + y(8)^n(ii,6)) - y(10))/tau(ii,10); 
                     ((w(ii,22)*(k(ii,22)^n(ii,22)+1)*y(7)^n(ii,22)*ymax)/(k(ii,22)^n(ii,22) + y(7)^n(ii,22)) - y(11))/tau(ii,11); 
                     (w(ii,25)*(k(ii,25)^n(ii,25)+1)*y(11)^n(ii,25)/(k(ii,25)^n(ii,25) + y(11)^n(ii,25)) + w(ii,34)*(k(ii,34)^n(ii,34)+1)*y(22)^n(ii,34)/(k(ii,34)^n(ii,34) + y(22)^n(ii,34)) - (w(ii,25)*(k(ii,25)^n(ii,25)+1)*y(11)^n(ii,25)/(k(ii,25)^n(ii,25) + y(11)^n(ii,25)))*(w(ii,34)*(k(ii,34)^n(ii,34)+1)*y(22)^n(ii,34)/(k(ii,34)^n(ii,34) + y(22)^n(ii,34)))*ymax - y(12))/tau(ii,12); 
                     (w(ii,9)*(k(ii,9)^n(ii,9)+1)*y(10)^n(ii,9)/(k(ii,9)^n(ii,9) + y(10)^n(ii,9)) + w(ii,11)*(k(ii,11)^n(ii,11)+1)*y(14)^n(ii,11)/(k(ii,11)^n(ii,11) + y(14)^n(ii,11)) - ((w(ii,9)*(k(ii,9)^n(ii,9)+1)*y(10)^n(ii,9)/(k(ii,9)^n(ii,9) + y(10)^n(ii,9)))*(w(ii,11)*(k(ii,11)^n(ii,11)+1)*y(14)^n(ii,11)/(k(ii,11)^n(ii,11) + y(14)^n(ii,11))))*ymax - y(13))/tau(ii,13); 
                     ((w(ii,8)*(k(ii,8)^n(ii,8)+1)*y(9)^n(ii,8)/(k(ii,8)^n(ii,8) + y(9)^n(ii,8)))*ymax - y(14))/tau(ii,14); 
                     ((w(ii,10)*(k(ii,10)^n(ii,10)+1)*y(14)^n(ii,10)/(k(ii,10)^n(ii,10) + y(14)^n(ii,10)))*ymax - y(15))/tau(ii,15); 
                     (w(ii,23)*(k(ii,23)^n(ii,23)+1)*y(5)^n(ii,23)/(k(ii,23)^n(ii,23) + y(5)^n(ii,23)) + w(ii,24)*(k(ii,24)^n(ii,24)+1)*y(4)^n(ii,24)/(k(ii,24)^n(ii,24) + y(4)^n(ii,24)) - ((w(ii,23)*(k(ii,23)^n(ii,23)+1)*y(5)^n(ii,23)/(k(ii,23)^n(ii,23) + y(5)^n(ii,23)))*(w(ii,24)*(k(ii,24)^n(ii,24)+1)*y(4)^n(ii,24)/(k(ii,24)^n(ii,24) + y(4)^n(ii,24))))*ymax - y(16))/tau(ii,16); 
                     ((w(ii,26)*(k(ii,26)^n(ii,26)+1)*y(16)^n(ii,26)/(k(ii,26)^n(ii,26) + y(16)^n(ii,26)))*ymax - y(17))/tau(ii,17); 
                     (w(ii,30)*(k(ii,30)^n(ii,30)+1)*y(12)^n(ii,30)/(k(ii,30)^n(ii,30) + y(12)^n(ii,30)) +  w(ii,29)*(k(ii,29)^n(ii,29)+1)*y(26)^n(ii,29)/(k(ii,29)^n(ii,29) + y(26)^n(ii,29)) - ((w(ii,30)*(k(ii,30)^n(ii,30)+1)*y(12)^n(ii,30)/(k(ii,30)^n(ii,30) + y(12)^n(ii,30)))*(w(ii,29)*(k(ii,29)^n(ii,29)+1)*y(26)^n(ii,29)/(k(ii,29)^n(ii,29) + y(26)^n(ii,29))))*ymax - y(18))/tau(ii,18); 
                     ((w(ii,7)*(k(ii,7)^n(ii,7)+1)*y(9)^n(ii,7)/(k(ii,7)^n(ii,7) + y(9)^n(ii,7))*w(ii,7)*(k(ii,7)^n(ii,7)+1)*y(13)^n(ii,7)/(k(ii,7)^n(ii,7) + y(13)^n(ii,7))) + w(ii,13)*(k(ii,13)^n(ii,13)+1)*y(15)^n(ii,13)/(k(ii,13)^n(ii,13) + y(15)^n(ii,13)) - ((w(ii,7)*(k(ii,7)^n(ii,7)+1)*y(9)^n(ii,7)/(k(ii,7)^n(ii,7) + y(9)^n(ii,7))*w(ii,7)*(k(ii,7)^n(ii,7)+1)*y(13)^n(ii,7)/(k(ii,7)^n(ii,7) + y(13)^n(ii,7))))*(w(ii,13)*(k(ii,13)^n(ii,13)+1)*y(15)^n(ii,13)/(k(ii,13)^n(ii,13) + y(15)^n(ii,13)))*ymax - y(19))/tau(ii,19); 
                     (w(ii,33)*(k(ii,33)^n(ii,33)+1)*y(22)^n(ii,33)/(k(ii,33)^n(ii,33) + y(22)^n(ii,33)) + w(ii,40)*(k(ii,40)^n(ii,40)+1)*y(29)^n(ii,40)/(k(ii,40)^n(ii,40) + y(29)^n(ii,40)) - ((w(ii,33)*(k(ii,33)^n(ii,33)+1)*y(29)^n(ii,33)/(k(ii,33)^n(ii,33) + y(29)^n(ii,33)))*(w(ii,40)*(k(ii,40)^n(ii,40)+1)*y(29)^n(ii,40)/(k(ii,40)^n(ii,40) + y(29)^n(ii,40))))*ymax - y(20))/tau(ii,20);
                     (w(ii,35)*(k(ii,35)^n(ii,35)+1)*y(12)^n(ii,35)/(k(ii,35)^n(ii,35) + y(12)^n(ii,35))*(w(ii,35)*(k(ii,35)^n(ii,35)+1)*y(20)^n(ii,35)/(k(ii,35)^n(ii,35) + y(20)^n(ii,35)))*ymax - y(21))/tau(ii,21); 
                     (w(ii,27)*(k(ii,27)^n(ii,27)+1)*y(17)^n(ii,27)/(k(ii,27)^n(ii,27) + y(17)^n(ii,27))*ymax - y(22))/tau(ii,22); 
                     ((w(ii,14)*(k(ii,14)^n(ii,14)+1)*y(19)^n(ii,14)/(k(ii,14)^n(ii,14) + y(19)^n(ii,14))) + (w(ii,31)*(k(ii,31)^n(ii,31)+1)*y(18)^n(ii,31)/(k(ii,31)^n(ii,31) + y(18)^n(ii,31))) - ((w(ii,31)*(k(ii,31)^n(ii,31)+1)*y(18)^n(ii,31)/(k(ii,31)^n(ii,31) + y(18)^n(ii,31)))*(w(ii,14)*(k(ii,14)^n(ii,14)+1)*y(19)^n(ii,14)/(k(ii,14)^n(ii,14) + y(19)^n(ii,14))))*ymax - y(23))/tau(ii,23); 
                     ((w(ii,12)*(k(ii,12)^n(ii,12)+1)*y(18)^n(ii,12)/(k(ii,12)^n(ii,12) + y(18)^n(ii,12))) + (w(ii,15)*(k(ii,15)^n(ii,15)+1)*y(19)^n(ii,15)/(k(ii,15)^n(ii,15) + y(19)^n(ii,15))) - ((w(ii,15)*(k(ii,15)^n(ii,15)+1)*y(19)^n(ii,15)/(k(ii,15)^n(ii,15) + y(19)^n(ii,15)))*(w(ii,12)*(k(ii,12)^n(ii,12)+1)*y(18)^n(ii,12)/(k(ii,12)^n(ii,12) + y(18)^n(ii,12))))*ymax - y(24))/tau(ii,24); 
                     ((w(ii,18)*(k(ii,18)^n(ii,18)+1)*y(19)^n(ii,18)/(k(ii,18)^n(ii,18) + y(19)^n(ii,18))) + (w(ii,32)*(k(ii,32)^n(ii,32)+1)*y(18)^n(ii,32)/(k(ii,32)^n(ii,32) + y(18)^n(ii,32))) - ((w(ii,18)*(k(ii,18)^n(ii,18)+1)*y(19)^n(ii,18)/(k(ii,18)^n(ii,18) + y(19)^n(ii,18)))*(w(ii,32)*(k(ii,32)^n(ii,32)+1)*y(18)^n(ii,32)/(k(ii,32)^n(ii,32) + y(18)^n(ii,32))))*ymax - y(25))/tau(ii,25); 
                     (w(ii,28)*(k(ii,28)^n(ii,28)+1)*y(4)^n(ii,28)/(k(ii,28)^n(ii,28) + y(4)^n(ii,28))*ymax - y(26))/tau(ii,26); 
                     (w(ii,17)*(k(ii,17)^n(ii,17)+1)*y(6)^n(ii,17)/(k(ii,17)^n(ii,17) + y(6)^n(ii,17))*ymax - y(27))/tau(ii,27);
                     (w(ii,38)*(k(ii,38)^n(ii,38)+1)*y(29)^n(ii,38)/(k(ii,38)^n(ii,38) + y(29)^n(ii,38))*ymax - y(28))/tau(ii,28); 
                     ((1-w(ii,36)*(k(ii,36)^n(ii,36)+1)*y(20)^n(ii,36)/(k(ii,36)^n(ii,36) + y(20)^n(ii,36))) + (w(ii,37)*(k(ii,37)^n(ii,37)+1)*y(26)^n(ii,37)/(k(ii,37)^n(ii,37) + y(26)^n(ii,37))) - (w(ii,37)*(k(ii,37)^n(ii,37)+1)*y(26)^n(ii,37)/(k(ii,37)^n(ii,37) + y(26)^n(ii,37)))*((1-w(ii,36)*(k(ii,36)^n(ii,36)+1)*y(20)^n(ii,36)/(k(ii,36)^n(ii,36) + y(20)^n(ii,36))))*ymax - y(29))/tau(ii,29); 
                     (w(ii,39)*(k(ii,39)^n(ii,39)+1)*y(28)^n(ii,39)/(k(ii,39)^n(ii,39) + y(28)^n(ii,39))*ymax - y(30))/tau(ii,30);
                     ];
   
    
    % solve using numerical ODE solver 
    [t,sol] = ode23s(diffEq,[firstTime:lastTime],y0,odeOpts); % [0:48] timesteps instead of [0 48]
    % interpolate solution to specified timesteps
    ROS = interp1(t,sol(:,13),48);
    IL6 = interp1(t,sol(:,23),48);
    TNFa = interp1(t,sol(:,24),48);
    IL1b = interp1(t,sol(:,25),48);
    VEGFa = interp1(t,sol(:,27),48);
    ROSec = interp1(t,sol(:,12),48);
    NO = interp1(t,sol(:,20),48);
    eNOS = interp1(t,sol(:,22),48);

    % assign solution
    Y(ii,:) = [ROS' , IL6', TNFa', IL1b',VEGFa',ROSec',NO',eNOS'];

end



