function azim = VectMean_arctan(dx, dy)

% VectMean_arctan.m

% Copyright C 2004  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.

% VectMean_arctan returns arctangent (in degrees) of angle defined by dx/dy
%   arctangent is set up assuming AZIMUTH data in which North=0, East=90
%   dx and dy are shifts in the X and Y direction
%   ===>>> (sum sines and sum cosines , respectively)

% process for 0, 90, 180, 270

ang = -999;
if (dx == 0 & dy > 0)
    ang = 0;
end    
if (dx == 0 & dy < 0)
    ang = 180;
end        
if (dx > 0 & dy == 0)
    ang = 90;
end        
if (dx < 0 & dy == 0)
    ang = 270; 
end
if  ang > -99
    azim = ang;
    return
end    

% process for non-axis directions 

ang = atan(abs(dx/dy))*57.3;
    if (dx > 0 & dy < 0)
        ang = 180 - ang;
    end    
    if (dx < 0 & dy > 0)
        ang = 360 - ang;
    end    
    if (dx < 0 & dy < 0)
        ang = ang + 180;
    end
azim = ang;
