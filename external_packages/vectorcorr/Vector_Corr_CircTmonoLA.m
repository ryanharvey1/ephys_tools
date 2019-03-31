function RtrnCode = Vector_Corr_CircTmonoLA(AzVect1, AzVect2, testalfa, ...
                                            fidO)

% Vector_Corr_CircTmonoLA.m

% Copyright C 2005  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.
       
% circular-circular correlation  -  PIhatN                                  
% Distribution-free rank-correlation method
% Uses large-sample approx. for distribution of test
% Fisher, 1993, p. 148 - 149
% Mardia and Jupp, 2000, p. 252

% Input variables:
%   AzVect1: vector of azimuths, variable 1
%   AzVect2: vector of azimuths, variable 2
%   testalfa: significance level for tests of significant correlation
%          only tabled for 0.10, 0.05, 0.025, 0.01, 0.005
%   fidO: output unit for external file (write if fidO>0)
% Output variable:
%   RtrnCode: return code

% Functions called from this module:
%   Vector_Corr_CirRankPiTest.m
                                  
RtrnCode = 0;

% find nearest tabled value to specified alfa
alfalist = [0.10, 0.05, 0.025, 0.01, 0.005];
tblalfa = alfalist(1); diff = abs(testalfa - alfalist(1));
for iii = 2:5
    dd = abs(testalfa - alfalist(iii));
    if dd < diff
        diff = dd; tblalfa = alfalist(iii);
    end
end

NTot = length(AzVect1);

     % get ranks of two variables
          
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

if NTot >= 8
  U = (NTot - 1)*abs(PIhatN);
  alfa = tblalfa;
  cutoff = Vector_Corr_CirRankPiTest(NTot, alfa);
end

     % output results
 
pt20 = sprintf('   Large-sample approximation to test distribution\n');       
pt21 = sprintf('     Estimated correlation (PIhatN) = %7.3g\n', PIhatN);
pt22B= sprintf('     No test of significance: invalid alfa = %.5f\n',alfa);
pt22 = sprintf('     Test statistic = %9.5g\n', U);
pt22alf = sprintf('     Significance level: alfa = %.3f\n', alfa);  
pt23 = sprintf('     Test criterion (cutoff) = %9.5g\n', cutoff);
if U > cutoff
    pt24 = sprintf('     Reject hypothesis of no association\n');
else
    pt24 = sprintf('     Cannot reject hypothesis of no association\n');
end
ptskip = sprintf('\n');

disp(pt20); disp(pt21); 
if NTot >= 8 
    if cutoff > -990
       disp(pt22); disp(pt22alf); disp(pt23); disp(pt24);
    else
       disp(pt22B);
    end
end
if NTot < 25
    pt26 = sprintf('     NOTE: Test is approximate due to small sample\n');
    disp(pt26)
end
disp(ptskip);
if fidO > 0 
    fprintf(fidO, pt20); fprintf(fidO, pt21); 
    if NTot >= 8 
        if cutoff > -990
          fprintf(fidO, pt22); fprintf(fidO, pt22alf);
          fprintf(fidO, pt23); fprintf(fidO, pt24);
        else
          fprintf(fidO, pt22B);
        end
    end
    if NTot < 25 
        fprintf(fidO, pt26);
    end
    fprintf(fidO, ptskip);
end

clear xxx Indx RankAz1 RankAz2 rankcos1 ranksin1 rankcos2 ranksin2