function RtrnCode = Vector_Corr_CircCirc1(AzVect1, AzVect2, alfa, fidO, ...
                                          NTrials)

% Vector_Corr_CircCirc1.m

% Copyright C 2005  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.

% circular-circular  -  R_squared
% Treats both azimuth variables as bivariate (sin and cos components)
% Then gets multiple canonical correlations, and calcs. rsquared = sum
% Embedded method with multiple correlation on 4 variables
% Large-sample approx or resampling method
% Mardia and Jupp, 2000, p. 248-249
% Resampling: Fisher, 1993, p, 214 - 218

% Input variables:
%   AzVect1: vector of azimuths, variable 1
%   AzVect2: vector of azimuths, variable 2
%   alfa: significance level for tests of significant correlation
%   fidO: output unit for external file (write if fidO>0)
%   NTrials: number of resampling trials
% Output variable:
%   RtrnCode: return code

% Functions used:
%   Vector_Corr_CircCirc1LA.m
%   Vector_Corr_CircCirc1R.m
                                  
RtrnCode = 0;
                       
pt11 = sprintf(['   Sum of embedded (4-variate) multiple canonical', ...
                     ' correlations\n', ...
                '   Ref.: Mardia and Jupp, 2000, p. 248 - 249\n']);    
disp(pt11)
if fidO > 0 
   fprintf(fidO, pt11);
end

NTot = length(AzVect1);

if NTot > 14
    
    % process assuming large-sample approximation -
    % overlapping resampling for N = 15 - 30
    
    RtrnCode = Vector_Corr_CircCirc1LA(AzVect1, AzVect2, alfa, fidO);
    
end

if NTot < 31
    
    % process assuming resampling (permutation test) -
    % overlapping large-sample approx. for N = 15 - 30
    
    RtrnCode = Vector_Corr_CircCirc1R(AzVect1, AzVect2, alfa, fidO, ...
                                    NTrials);
    
end
 
ptlin = sprintf('     -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -\n');
ptskip = sprintf('\n');
disp(ptlin); disp(ptskip)
if fidO > 0 
    fprintf(fidO, ptlin); fprintf(fidO, ptskip);
end

