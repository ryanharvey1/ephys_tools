function [RtrnCode, FigCntX] = Vector_Corr_Plots(type, Vect1, Vect2, ...
                   xlbl, ylbl, DataTtl, FigCnt)    

%Vector_Corr_Plots.m

% Copyright C 2004  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.

% Generate crossplots of pair of variables for correlation
% For azimuthal variables, doubles range of angles to [0,720] and
%           repeats the angles + 360 degrees

% Input variables:
%   type: type of crossplot
%   1 = azimuth-azimuth       2 = linear-azimuth   
%   Vect1: Data for Y axis (azimuths for types 1 and 2) - degrees
%   Vect2: Data for X axis (linear for type=2; azimuth for type=1: degrees)
%   xlbl: label for X axis
%   ylbl: label for Y axis
%   DataTtl: Title for data or job
%   FigCnt: Counter to number figures as are generated
% Output variable:
%   RtrnCode: Error return code
%   FigCntX: Number of figures generated (for naming next figure)
   

RtrnCode = 0;
NTot = length(Vect1);
FigCnt = FigCnt + 1;
figure(FigCnt)

if type == 1
    
    % circular-circular plot
    
Xvect = zeros(4*NTot, 1);
Xvect(1:NTot,1) = Vect1;
Xvect(NTot+1:2*NTot,1) = Vect1 + 360;
Xvect(2*NTot+1:3*NTot,1) = Vect1;
Xvect(3*NTot+1:4*NTot,1) = Vect1 + 360;
Yvect = zeros(4*NTot, 1);
Yvect(1:NTot) = Vect2;  
Yvect(NTot+1:2*NTot) = Vect2 + 360; 
Yvect(2*NTot+1:3*NTot,1) = Vect2 + 360;
Yvect(3*NTot+1:4*NTot,1) = Vect2 ;
       
plot(Xvect, Yvect, 'x')
V = axis;
axis([0 720  0 720])
xlbl2 = ['Azimuth: ', xlbl];
ylbl2 = ['Azimuth: ', ylbl];
xlabel(xlbl2), ylabel(ylbl2)
line([0 720] , [360 360], 'LineStyle', '-.')
line([360 360] , [0 720], 'LineStyle', '-.')
pthdr2=sprintf('%s', DataTtl);  
title(pthdr2, 'FontSize', 12);

else
 
    % circular-linear plot

Xvect = zeros(2*NTot, 1);
Xvect(1:NTot,1) = Vect2;
Xvect(NTot+1:2*NTot,1) = Vect2;
Yvect = zeros(2*NTot, 1);
Yvect(1:NTot) = Vect1;  
Yvect(NTot+1:2*NTot) = Vect1 + 360; 
        
plot(Xvect, Yvect, 'x')
V = axis;
axis([V(1) V(2)  0 720])
ylbl2 = ['Azimuth: ', ylbl];
xlabel(xlbl), ylabel(ylbl2)
line([V(1) V(2)] , [360 360], 'LineStyle', '-.')
pthdr2=sprintf('%s', DataTtl);  
title(pthdr2, 'FontSize', 12);
    
end

FigCntX = FigCnt;
        
clear Xvect Yvect        
        