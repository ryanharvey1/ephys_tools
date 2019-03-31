function RtrnCode = Vector_Corr_CircCirc1R(AzVect1, AzVect2, alfa, ...
                                           fidO, NTrials)

% Vector_Corr_CircCirc1R.m

% Copyright C 2005  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.

% circular-circular  -  R_squared
% Treats both azimuth variables as bivariate (sin and cos components)
% Then gets multiple canonical correlations, and calcs. rsquared=sum
% Embedded method with multiple correlation on 4 variables
% Uses resampling to get test distribution
% Mardia and Jupp, 2000, p. 248-249
% Resampling: Fisher, 1993, p. 214 - 218

% Input variables:
%   AzVect1: vector of azimuths, variable 1
%   AzVect2: vector of azimuths, variable 2
%   alfa: significance level for tests of significant correlation
%   fidO: output unit for external file (write if fidO>0)
%   NTrials: number of resampling trials
% Output variable:
%   RtrnCode: return code
                               
RtrnCode = 0;

NTot = length(AzVect1);

% get basic calculations, including estimate of coefficient

     % get single correlations to calculate multiple canonical correlation
     
DArray = zeros(NTot, 4);
AzV1 = AzVect1/57.3;   AzV2 = AzVect2/57.3;
DArray(:,1) = cos(AzV1);
DArray(:,2) = sin(AzV1);
DArray(:,3) = cos(AzV2);
DArray(:,4) = sin(AzV2);

% alternative way to calculate r-squared - Mardia, Kent, & Bibby
%CovArray = cov(DArray)
%S11 = CovArray(1:2, 1:2); %S22 = CovArray(3:4, 3:4)
%S12 = CovArray(1:2, 3:4); %S21 = transpose(S12)
%prod = inv(S11)*S12*inv(S22)*S21
%rSq = trace(prod)

RArray = corrcoef(DArray);
r1 = RArray(1,2);
r2 = RArray(3,4);
rcc = RArray(1,3);
rcs = RArray(1,4);
rsc = RArray(2,3);
rss = RArray(2,4);

R_Sq = rcc^2 + rcs^2 + rsc^2 + rss^2 + 2*(rcc*rss + rcs*rsc)*r1*r2;
R_Sq = R_Sq - 2*(rcc*rcs + rsc*rss)*r2 - 2*(rcc*rsc + rcs*rss)*r1;
R_Sq = R_Sq/((1 - r1^2)*(1 - r2^2));

% set up and loop over NTrials permutation trials for resampling

SaveStat = zeros(1, NTrials);
AzV3 = zeros(1, NTot);

for iii = 1:NTrials
        
    % get random permutation of linear variables
    
    R = randperm(NTot);
    
    for jjj = 1:NTot    
       AzV3(jjj) = AzV2(R(jjj));   
    end
     
     % calculate stat for measure-of-association coefficient
     
     DArray(:,3) = transpose(cos(AzV3));
     DArray(:,4) = transpose(sin(AzV3));

     RArray = corrcoef(DArray);
     r1 = RArray(1,2);
     r2 = RArray(3,4);
     rcc = RArray(1,3);
     rcs = RArray(1,4);
     rsc = RArray(2,3);
     rss = RArray(2,4);

     R2 = rcc^2 + rcs^2 + rsc^2 + rss^2 + 2*(rcc*rss + rcs*rsc)*r1*r2;
     R2 = R2 - 2*(rcc*rcs + rsc*rss)*r2 - 2*(rcc*rsc + rcs*rss)*r1;
     R2 = R2/((1 - r1^2)*(1 - r2^2));
     SaveStat(iii) = R2;
     
end     
     
% sort calculated statistics and get cutoffs from generated distrib.

SortStat = sort(SaveStat);

L1 = fix(NTrials*(1 - alfa) + 0.5);
cutoff = SortStat(L1);   

     % output results
     
pt20 = sprintf(['   Test distribution based on resampling\n', ...
                '   Ref.: Fisher, 1993, p. 214 - 218\n']); 
pt21 = sprintf('     Estimated sum correlations (Rsquared) = %7.3g\n',R_Sq);
pt22a = sprintf('     Significance level: alfa = %.3f\n', alfa);  
pt23 = sprintf('     Test criterion (cutoff) = %7.3g\n', cutoff);
if R_Sq > cutoff
    pt24 = sprintf('     Reject hypothesis of no association\n');
else
    pt24 = sprintf('     Cannot reject hypothesis of no association\n');
end

ptskip = sprintf('\n');

disp(pt20); disp(pt21); disp(pt22a); disp(pt23); disp(pt24); 
disp(ptskip)
if fidO > 0 
    fprintf(fidO, pt20); fprintf(fidO, pt21); fprintf(fidO, pt22a);
    fprintf(fidO, pt23); fprintf(fidO, ptskip);
end

clear RArray DArray AzV1 AzV2 AzV3 SaveStat SortStat R iii jjj 