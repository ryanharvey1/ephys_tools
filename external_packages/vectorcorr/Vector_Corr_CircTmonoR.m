function RtrnCode = Vector_Corr_CircTmonoR(AzVect1, AzVect2, alfa, ...
                                           fidO, NTrials)

% Vector_Corr_CircTmonoR.m

% Copyright C 2005  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.
       
% circular-circular correlation  -  PIhatN                                  
% Distribution-free rank-correlation method
% Uses resampling for test distribution
% Fisher, 1993, p. 148 - 149
% Mardia and Jupp, 2000, p. 252
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

     % get ranks of two azimuthal variables
          
[xxx, Indx] = sort(AzVect1);
jjj = 0;
for iii = 1:NTot
    jjj = jjj + 1;
    RankAz1(Indx(iii)) = jjj;
end

[xxx, Indx] = sort(AzVect2);
jjj = 0;
for iii = 1:NTot
    jjj = jjj + 1;
    RankAz2(Indx(iii)) = jjj;
end

     % calculate measure of association and test significance
     
RankAz1 = 2*pi*RankAz1/NTot;     RankAz2 = 2*pi*RankAz2/NTot;
rankcos1 = cos(RankAz1);         ranksin1 = sin(RankAz1);
rankcos2 = cos(RankAz2);         ranksin2 = sin(RankAz2);

A = rankcos1*transpose(rankcos2);
B = ranksin1*transpose(ranksin2);
C = rankcos1*transpose(ranksin2);
D = ranksin1*transpose(rankcos2);

PIhatN = 4*(A*B - C*D)/NTot^2;

% set up and loop over permutation trials

SaveStat = zeros(1, NTrials);
AzVect3 = zeros(1, NTot);

for iii = 1:NTrials
        
    % get random permutation of linear variables
    
    R = randperm(NTot);
    
    for jjj = 1:NTot    
       AzVect3(jjj) = AzVect2(R(jjj));   
    end
     
     % calculate stat for measure-of-association coefficient

    [xxx, Indx] = sort(AzVect3);
    ijjj = 0;
    for iiii = 1:NTot
        ijjj = ijjj + 1;
        RankAz2(Indx(iiii)) = ijjj;
    end

    RankAz2 = 2*pi*RankAz2/NTot;
    rankcos2 = cos(RankAz2);         ranksin2 = sin(RankAz2);

    A = rankcos1*transpose(rankcos2);
    B = ranksin1*transpose(ranksin2);
    C = rankcos1*transpose(ranksin2);
    D = ranksin1*transpose(rankcos2);

    SaveStat(iii) = A*B - C*D;
    
end

% sort calculated statistics and get cutoffs from generated distrib.

SortStat = sort(SaveStat);

L2 = fix(NTrials*(1 - alfa/2) + 0.5);
cutoff2 = 4*SortStat(L2)/NTot^2;
L1 = NTrials - L2;
cutoff1 = 4*SortStat(L1)/NTot^2;

     % output results
 
pt20 = sprintf(['   Test distribution based on resampling\n', ...
                '   Ref.: Fisher, 1993, p. 214 - 218\n']);       
pt21 = sprintf('     Estimated correlation (PIhatN) = %7.3g\n', PIhatN);
pt22alf = sprintf('     Significance level: alfa = %.3f\n', alfa);  
pt23 = sprintf('     Test criteria (cutoffs) = %7.3g,%7.3g\n', ...
                                cutoff1, cutoff2);
if PIhatN < cutoff1 | PIhatN > cutoff2
    pt24 = sprintf('     Reject hypothesis of no association\n');
else
    pt24 = sprintf('     Cannot reject hypothesis of no association\n');
end
ptskip = sprintf('\n');

disp(pt20); disp(pt21); 
disp(pt22alf); disp(pt23); disp(pt24);
disp(ptskip)
if fidO > 0 
    fprintf(fidO, pt20); fprintf(fidO, pt21); 
    fprintf(fidO, pt22alf);
    fprintf(fidO, pt23); fprintf(fidO, pt24);
    fprintf(fidO, ptskip);
end

clear xxx Indx RankAz1 RankAz2 rankcos1 ranksin1 rankcos2 ranksin2
clear iii jjj iiii ijjj AzVect3 SaveStat SortStat R
