function Gp = step_function(t, glu_sampled)

global number_ctrl time_ctrl density diameter GC_conc GC_time GC_LB GC_UB time_lee glu_UB


% poly1/linear fit before 12 weeks
% l(t) = 0.021*t + 10.33

% poly2 fit before 12 weeks
% f(t) = -2.994e-05*t^2 + 0.08638*t - 18.4


% %Run the glucose minimal model first and get max and min values, then run it again with normalized lbm input
% % if (t>=336 && t < time_finch(1))
% %         Gp = 0.021*t + 10.33;
% %         %Gp = -2.994e-05*t^2 + 0.08638*t - 18.4;
% % else
% syms t
% Gp(t) = piecewise(t > GC_time(1) & t < GC_time(2), glu_sampled(1,1),...
%                t >= GC_time(2) & t < GC_time(3),glu_sampled(2,1), ...
%                t >= GC_time(3) & t < GC_time(4),glu_sampled(3,1), ...
%                t >= GC_time(4) & t < GC_time(5),glu_sampled(4,1), ...
%                t >= GC_time(5) & t < GC_time(6),glu_sampled(5,1), ...
%                t >= GC_time(6) & t < GC_time(7),glu_sampled(6,1), ...
%                t >= GC_time(7) & t < GC_time(8),glu_sampled(7,1), ...
%                t >= GC_time(8) & t < GC_time(9),glu_sampled(8,1), ...
%                t >= GC_time(9) & t < GC_time(10),glu_sampled(9,1), ...
%                t >= GC_time(10) & t < GC_time(11),glu_sampled(10,1), ...
%                glu_sampled(11,1));


if (t > GC_time(1) && t < GC_time(2))
        Gp = glu_sampled(1,1);
elseif (t >= GC_time(2)) && (t < GC_time(3))
        Gp = glu_sampled(2,1);  

elseif (t >= GC_time(3)) && (t < GC_time(4))
        Gp = glu_sampled(3,1);  

elseif (t >= GC_time(4)) && (t < GC_time(5))
        Gp = glu_sampled(4,1); 

elseif (t >= GC_time(5)) && (t < GC_time(6))
        Gp = glu_sampled(5,1);
        
elseif (t >= GC_time(6)) && (t < GC_time(7))
        Gp = glu_sampled(6,1);

elseif (t >= GC_time(7)) && (t < GC_time(8))
        Gp = glu_sampled(7,1);

elseif (t >= GC_time(8)) && (t < GC_time(9))
        Gp = glu_sampled(8,1);

elseif (t >= GC_time(9)) && (t < GC_time(10))
        Gp = glu_sampled(9,1);

elseif (t >= GC_time(10)) && (t < GC_time(11))
        Gp = glu_sampled(10,1);

elseif (t >= GC_time(11))
        Gp = glu_sampled(11,1);
end

% if (t > 840 && t < time_finch(1))
%         Gp = 41.2;
% elseif (t >= time_finch(1)) && (t < time_finch(2))
%         Gp = glu_sampled(1,1);  
% 
% elseif (t >= time_finch(2)) && (t < time_finch(3))
%         Gp = glu_sampled(2,1);  
% 
% elseif (t >= time_finch(3)) && (t < time_finch(4))
%         Gp = glu_sampled(3,1); 
% 
% elseif (t >= time_finch(4)) && (t < time_finch(5))
%         Gp = glu_sampled(4,1); 
%         
% elseif (t >= time_finch(5))
%         Gp = glu_sampled(5,1);
% 
% 
% end
