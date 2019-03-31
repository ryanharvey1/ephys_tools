function RtrnCode = Vector_Corr_CircTmono(AzVect1, AzVect2, testalfa, ...
                                          fidO, NTrials)

% Vector_Corr_CircTmono.m

% Copyright C 2005  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.
       
% circular-circular correlation  -  PIhatN                                  
% Distribution-free rank-correlation method
% Uses large-sample approx. or resampling
% Fisher, 1993, p. 148 - 149
% Mardia and Jupp, 2000, p. 252
% Resampling: Fisher, 1993, p. 214 - 218

% Input variables:
%   AzVect1: vector of azimuths, variable 1
%   AzVect2: vector of azimuths, variable 2
%   testalfa: significance level for tests of significant correlation
%          only tabled for 0.10, 0.05, 0.025, 0.01, 0.005
%   fidO: output unit for external file (write if fidO>0)
%   NTrials: number of resampling trials
% Output variable:
%   RtrnCode: return code

% Functions called from this module:
%   Vector_Corr_CircTmonoLA.m
%   Vector_Corr_CircTmonoR.m
                                  
RtrnCode = 0;
                             
pt11 = sprintf(['   T-monotonic rank-correlation association\n', ...
                '   Note: Test is sensitive to ties in ', ...
                    'azimuth data\n', ...
                '   Refs.: Fisher, 1993, p. 148 - 149\n', ...   
                '          Mardia and Jupp, 2000, p. 252\n']);    
disp(pt11)
if fidO > 0 
   fprintf(fidO, pt11);
end

NTot = length(AzVect1);

if NTot > 7
    
    % process assuming large-sample approximation -
    % overlapping resampling for N = 8 - 30
    
    RtrnCode = Vector_Corr_CircTmonoLA(AzVect1, AzVect2, testalfa, fidO);
    
end

if NTot < 31
    
    % process assuming resampling (permutation test) -
    % overlapping large-sample approx. for N = 8 - 30
    
    RtrnCode = Vector_Corr_CircTmonoR(AzVect1, AzVect2, testalfa, fidO,...
                                      NTrials);
    
end

ptlin = sprintf('     -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - \n');
ptskip = sprintf('\n');

disp(ptlin); disp(ptskip)
if fidO > 0 
    fprintf(fidO, ptlin); fprintf(fidO, ptskip);
end
