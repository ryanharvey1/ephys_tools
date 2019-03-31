function RtrnCode = Vector_Corr_CAssoc1(AzVect, LinVect, testalfa, ...
                                        fidO, NTrials)

% Vector_Corr_CAssoc1.m

% Copyright C 2005  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.

% Linear-vector association   -   Un, Dn                               
% Distribution-free rank method - watch out: sensitive to tied data
% Uses large-sample approx or resampling
% Mardia and Jupp, 2000, p. 246-248
% Fisher, 1993, p. 140-141
% Resampling: Fisher, 1993, p. 214 - 218

% Input variables:
%   AzVect: vector of azimuths
%   LinVect: vector of linear (X) variable 
%   testalfa: significance level for tests of significant correlation
%        tabled values are 0.10, 0.05, 0.01 for N < 100
%   fidO: output unit for external file (write if fidO>0)
%   NTrials: number of trials for resampling
% Output variable:
%   RtrnCode: return code

% Function called from this module:
%   Vector_Corr_CAssoc1LA.m
%   Vector_Corr_CAssoc1R.m
       
RtrnCode = 0;
                   
pt11 = sprintf(['   C-Association\n', ...
                '   Distribution-free rank correlation\n', ...
                '   Note: Test is sensitive to ties in ', ...
                    'X or Theta data\n', ...
                '   Refs.: Mardia and Jupp, 2000, p. 246-248\n',...
                '          Fisher, 1993, p. 140-141\n']);    
disp(pt11)
if fidO > 0 
   fprintf(fidO, pt11);
end

NTot = length(AzVect);

if NTot > 5
    
    % process assuming large-sample approximation -
    % overlapping resampling for N = 6 - 30
    
    RtrnCode = Vector_Corr_CAssoc1LA(AzVect, LinVect, testalfa, fidO);
    
end

if NTot < 31
    
    % process assuming resampling (permutation test) -
    % overlapping large-sample approx. for N = 6 - 30
    
    RtrnCode = Vector_Corr_CAssoc1R(AzVect, LinVect, testalfa, fidO, ...
                                    NTrials);
    
end
 
ptlin = sprintf('     -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -\n');
ptskip = sprintf('\n');
disp(ptlin); disp(ptskip)
if fidO > 0 
    fprintf(fidO, ptlin); fprintf(fidO, ptskip);
end
