%DEMDOAT  Auxiliary file used by DEMDOA
%         displays text and graphics describing the DOA problem.

echo off

% A. Swami April 15, 1993, Nov 15, 1994. 
% Copyright (c) 1991-2001 by United Signals & Systems, Inc.
%       $Revision: 1.7 $

%     RESTRICTED RIGHTS LEGEND
% Use, duplication, or disclosure by the Government is subject to
% restrictions as set forth in subparagraph (c) (1) (ii) of the 
% Rights in Technical Data and Computer Software clause of DFARS
% 252.227-7013. 
% Manufacturer: United Signals & Systems, Inc., P.O. Box 2374, 
% Culver City, California 90231. 
%
%  This material may be reproduced by or for the U.S. Government pursuant 
%  to the copyright license under the clause at DFARS 252.227-7013. 


clf

x = 0:10:70;   y = x *0; 
x00 = 35; 
plot(x,y,'w:', x,y,'g+'),                         % sensors
hold on

x = x(1:5); 

th = -(-25*pi/180);   y0 = x00 * tan(th); 
y = x*tan(th)+y0;  plot(x, y, 'r')                 % wavefronts 
x0 = -y(2)*cot(th+pi/2) + x(2); 
plot([x(2),x0], [y(2),0], 'y');  % normals 

th = -(-15*pi/180);   y0 = x00 * tan(th); 
y = x*tan(th)+y0; plot(x, y, 'r')                 % wavefronts 
x0 = -y(2)*cot(th+pi/2) + x(2); 
plot([x(2),x0], [y(2),0], 'y');  % normals 


th = -(30*pi/180);   y0 = x00 * tan(th); 
y = x*tan(th)+y0; 
x0 = -y(2)*cot(th+pi/2) + x(2); 
plot([x(2),x0], [y(2),0], 'b');  % normals 

th = -(-30*pi/180);   y0 = x00 * tan(th); 
y = x*tan(th)+y0; 
x0 = -y(2)*cot(th+pi/2) + x(2); 
plot([x(2),x0], [y(2),0], 'b');  % normals 


text('position',[-20,58,0], 'fontsize',12, 'color','g', ... 
     'fontweight','bold', ... 
     'string',' Source and sensor configuration for the DOA problem')
text('position',[-20,50,0],      'fontsize',10, ...
     'string', 'sensor locations denoted by +')
text('position',[-20,45,0], 'fontsize',10, ... 
     'string', 'source wave fronts depicted by red lines')
text('position',[-20,40,0], 'fontsize',10, ... 
     'string', 'direction of travel depicted by yellow lines')
text('position',[-20,35,0], 'fontsize',10, ... 
     'string', 'effective noise direction depicted in blue')


axis([-20 70 -15 55]) 
axis equal
axis off 
set (gcf, 'Name','HOSA: DOA Problem') 

hold off
return


