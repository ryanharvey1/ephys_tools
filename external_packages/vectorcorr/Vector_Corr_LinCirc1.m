function RtrnCode = Vector_Corr_LinCirc1(AzVect, LinVect, alfa, fidO, ...
                                         NTrials)

% Vector_Corr_LinCirc1.m

% Copyright C 2005  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.

% Mean(X) = function(Theta) 
% embedded bivariate correlation method
% Treats azimuth as bivariate variable (based on sin and cos components), 
% and then calculates joint correlation of these with X 
% Uses large-sample approx or resampling
% Mardia and Jupp, 2000, p. 245 - 246
% Fisher, 1993, p. 145
% Resampling: Fisher, 1993, p. 214 - 218

% Input variables:
%   AzVect: vector of azimuths
%   LinVect: vector of linear (X) variable 
%   alfa: significance level for tests of significant correlation
%   fidO: output unit for external file (write if fidO>0)
%   NTrials: number of resampling trials
% Output variable:
%   RtrnCode: return code

% Functions called:
%   Vector_Corr_LinCirc1LA.m
%   Vector_Corr_LinCirc1R.m
                
                                  
RtrnCode = 0;
                            
pt11 = sprintf(['   Mean(X) = function(Theta)\n', ...
                '   Embedded (bivariate) correlation\n', ...
                '   Ref.: Mardia and Jupp, 2000, p. 245-246\n', ...
                '         Fisher, 1993, p. 145\n']);    
disp(pt11)
if fidO > 0 
   fprintf(fidO, pt11);
end

NTot = length(AzVect);

if NTot > 19
    
    % process assuming large-sample approximation - 
    % overlapping resampling for N = 20 - 35
    
   RtrnCode = Vector_Corr_LinCirc1LA(AzVect, LinVect, alfa, fidO);
    
end

if NTot < 36
    
    % process using resampling (permutation test) method 
    % overlapping large-sample for N = 20 - 35
    
   RtrnCode = Vector_Corr_LinCirc1R(AzVect, LinVect, alfa, fidO, NTrials);
    
end

ptlin = sprintf('    -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -\n');
ptskip = sprintf('\n');

disp(ptlin); disp(ptskip)
if fidO > 0 
    fprintf(fidO, ptlin); fprintf(fidO, ptskip);
end
