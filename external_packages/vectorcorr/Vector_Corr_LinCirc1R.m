function RtrnCode = Vector_Corr_LinCirc1R(AzVect, LinVect, alfa, fidO, ...
                                          NTrials)

% Vector_Corr_LinCirc1R.m

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
% This uses resampling (permutation test) for the test distribution
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
                
                                  
RtrnCode = 0;

NTot = length(AzVect);

% get basic calculations, including estimate of coefficient
% use single correlations to calculate multiple correlation
     
DArray = zeros(NTot, 3);
DArray(:,1) = cos(AzVect/57.3);
DArray(:,2) = sin(AzVect/57.3);
DArray(:,3) = LinVect;

RArray = corrcoef(DArray);

rcs = RArray(1,2);
rxs = RArray(1,3);
rxc = RArray(2,3);

R_AL2 = (rxc^2 + rxs^2 - 2*rxc*rxs*rcs)/(1 - rcs^2); 

% set up and loop over number of permutation trials

SaveStat = zeros(1, NTrials);

for iii = 1: NTrials
    
    % get random permutation of linear variables
    
    R = randperm(NTot);
    
    % store randomized entries 
    
    for jjj = 1:NTot    
       DArray(jjj,3) = LinVect(R(jjj));   
    end
     
     % calculate stat for coefficient -  
     % get single correlations to calculate multiple correlation

     RArray = corrcoef(DArray);

     rcs = RArray(1,2);
     rxs = RArray(1,3);
     rxc = RArray(2,3);

     SaveStat(iii) = (rxc^2 + rxs^2 - 2*rxc*rxs*rcs)/(1 - rcs^2); 
       
end
       
% sort calculated statistics and get cutoffs from generated distribution

SortStat = sort(SaveStat);

L1 = fix(NTrials*(1 - alfa) + 0.5);
cutoff = SortStat(L1);

% output results
          
pt20A = sprintf(['   Test distribution based on resampling\n', ...
                '   Ref.: Fisher, 1993, p. 214 - 218\n']);     
pt21 = sprintf('     Estimated correlation (R2xa) = %7.3g\n', R_AL2);
pt22a = sprintf('     Significance level: alfa = %.3f\n', alfa);  
pt23 = sprintf('     Test criterion (cutoff) = %7.3g\n', cutoff);
if R_AL2 > cutoff
    pt24 = sprintf('     Reject hypothesis of no association\n');
else
    pt24 = sprintf('     Cannot reject hypothesis of no association\n');
end

ptskip = sprintf('\n');
disp(pt20A); 
disp(pt21); disp(pt22a); disp(pt23); disp(pt24); disp(ptskip)
if fidO > 0
    fprintf(fidO, pt20A); 
    fprintf(fidO, pt21); fprintf(fidO, pt22a);
    fprintf(fidO, pt23);
    fprintf(fidO, pt24); 
    fprintf(fidO, ptskip);
end

clear SaveStat SortStat R RArray DArray iii jjj


    
