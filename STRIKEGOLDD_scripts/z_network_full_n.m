
syms n3 n4 n5 n6 n7 n8 n9 n10 n11 n12 n13 n14 n15 n16 n17 n18 n19 n20 n21 n22 n23 n24 n25 n26 n27 n28 n29 n30 n31 n32 n33 n34 n35 n36 n37 n38 n39 n40...
    y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y12 y13  y14 y15 y16 y17 y18 y19 y20 y21 y22 y23 y24 y25 y26 y27 y28 y29 y30 ...
    
% states:
x    = [y1; y2; y3; y4; y5; y6; y7; y8; y9; y10; y11; y12; y13; y14; y15; y16; y17; y18; y19; y20; y21; y22; y23; y24; y25; y26; y27; y28; y29; y30];

% outputs:
h    = x;  

% no input:
u    = [];

% parameters:

p     = [n3; n4; n5; n6; n7; n8; n9; n10; n11; n12; n13; n14; n15; n16; n17; n18; n19; n20; n21; n22; n23; n24; n25; n26; n27; n28; n29; n30; n31; n32; n33; n34; n35; n36; n37; n38; n39; n40];

% known parameters
ymax = 1; w = 1; k = 0.5; t = 1;


% dynamic equations:    
f    = [(w*ymax - y1)/t;
        (w*ymax - y2)/t;
        ((w*(k^n4+1)*y1^n4*ymax)/(k^n4 + y1^n4) - y3)/t;
((w*(k^n19+1)*y27^n19*ymax)/(k^n19 + y27^n19)- y4)/t; 
((w*(k^n20+1)*y27^n20*ymax)/(k^n20 + y27^n20) - y5)/t;
((w*(k^n16+1)*y19^n16*ymax)/(k^n16 + y19^n16) - y6)/t; 
((w*(k^n21+1)*y3^n21*ymax)/(k^n21 + y3^n21) - y7)/t; 
((w*(k^n5+1)*y3^n5*ymax)/(k^n5 + y3^n5) - y8)/t; 
((w*(k^n3+1)*y2^n3*ymax)/(k^n3 + y2^n3) - y9)/t; 
((w*(k^n6+1)*y8^n6*ymax)/(k^n6 + y8^n6) - y10)/t; 
((w*(k^n22+1)*y7^n22*ymax)/(k^n22 + y7^n22) - y11)/t; 
(w*(k^n25+1)*y11^n25/(k^n25 + y11^n25) + w*(k^n34+1)*y22^n34/(k^n34 + y22^n34) - (w*(k^n25+1)*y11^n25/(k^n25 + y11^n25))*(w*(k^n34+1)*y22^n34/(k^n34 + y22^n34))*ymax - y12)/t; %OR 
(w*(k^n9+1)*y10^n9/(k^n9 + y10^n9) + w*(k^n11)+1*y14^n11/(k^n11 + y14^n11) - ((w*(k^n9+1)*y10^n9/(k^n9 + y10^n9))*(w*(k^n11+1)*y14^n11/(k^n11 + y14^n11)))*ymax - y13)/t; 
((w*(k^n8+1)*y9^n8/(k^n8 + y9^n8))*ymax - y14)/t; 
((w*(k^n10+1)*y14^n10/(k^n10 + y14^n10))*ymax - y15)/t; 
(w*(k^n23+1)*y5^n23/(k^n23 + y5^n23) + w*(k^n24+1)*y4^n24/(k^n24 + y4^n24) - ((w*(k^n23+1)*y5^n23/(k^n23 + y5^n23))*( w*(k^n24+1)*y4^n24/(k^n24 + y4^n24)))*ymax - y16)/t; 
((w*(k^n26+1)*y16^n26/(k^n26 + y16^n26))*ymax - y17)/t; 
(w*(k^n30+1)*y12^n30/(k^n30 + y12^n30) +  w*(k^n29+1)*y26^n29/(k^n29 + y26^n29) - ((w*(k^n30+1)*y12^n30/(k^n30 + y12^n30))*(w*(k^n29+1)*y26^n29/(k^n29 + y26^n29)))*ymax - y18)/t; 
((w*(k^n7+1)*y9^n7/(k^n7 + y9^n7)*w*(k^n7+1)*y13^n7/(k^n7 + y13^n7)) + w*(k^n13+1)*y15^n13/(k^n13 + y15^n13) - ((w*(k^n7+1)*y9^n7/(k^n7 + y9^n7)*w*(k^n7+1)*y13^n7/(k^n7 + y13^n7)))*(w*(k^n13+1)*y15^n13/(k^n13 + y15^n13))*ymax - y19)/t; 
(w*(k^n33+1)*y22^n33/(k^n33 + y22^n33) + w*(k^n40+1)*y29^n40/(k^n40 + y29^n40) - ((w*(k^n40+1)*y29^n40/(k^n40 + y29^n40))*(w*(k^n33+1)*y29^n33/(k^n33 + y29^n33)))*ymax - y20)/t;
(w*(k^n35+1)*y12^n35/(k^n35 + y12^n35)*(w*(k^n35+1)*y20^n35/(k^n35 + y20^n35))*ymax - y21)/t;
(w*(k^n27+1)*y17^n27/(k^n27 + y17^n27)*ymax - y22)/t; 
((w*(k^n14+1)*y19^n14/(k^n14 + y19^n14)) + (w*(k^n31+1)*y18^n31/(k^n31 + y18^n31)) - ((w*(k^n31+1)*y18^n31/(k^n31 + y18^n31))*(w*(k^n14+1)*y19^n14/(k^n14 + y19^n14)))*ymax - y23)/t; 
((w*(k^n12+1)*y18^n12/(k^n12 + y18^n12)) + (w*(k^n15+1)*y19^n15/(k^n15 + y19^n15)) - ((w*(k^n15+1)*y19^n15/(k^n15 + y19^n15))*(w*(k^n12+1)*y18^n12/(k^n12 + y18^n12)))*ymax - y24)/t; 
((w*(k^n18+1)*y19^n18/(k^n18 + y19^n18)) + (w*(k^n32+1)*y18^n32/(k^n32 + y18^n32)) - ((w*(k^n18+1)*y19^n18/(k^n18 + y19^n18))*(w*(k^n32+1)*y18^n32/(k^n32 + y18^n32)))*ymax - y25)/t; 
(w*(k^n28+1)*y4^n28/(k^n28 + y4^n28)*ymax - y26)/t; 
(w*(k^n17+1)*y6^n17/(k^n17 + y6^n17)*ymax - y27)/t;
(w*(k^n38+1)*y29^n38/(k^n38 + y29^n38)*ymax - y28)/t; 
((1-w*(k^n36+1)*y20^n36/(k^n36 + y20^n36)) + (w*(k^n37+1)*y26^n37/(k^n37 + y26^n37)) - (w*(k^n37+1)*y26^n37/(k^n37 + y26^n37))*((1-w*(k^n36+1)*y20^n36/(k^n36 + y20^n36)))*ymax - y29)/t; 
(w*(k^n39+1)*y28^n39/(k^n39 + y28^n39)*ymax - y30)/t;
]; 
    
% initial conditions:    
ics  = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

% which initial conditions are known:
known_ics = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1];

save('networkODE_n.mat','x','h','u','p','f','ics','known_ics');
