function RtrnCode = Vector_Corr_CAssoc1LA(AzVect, LinVect, testalfa, fidO)

% Vector_Corr_CAssoc1LA.m

% Copyright C 2005  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.

% Linear-vector association   -   Un, Dn                               
% Distribution-free rank method - watch out: sensitive to tied data
% Uses large-sample approximation for test distribution
% Mardia and Jupp, 2000, p. 246-248
% Fisher, 1993, p. 140-141

% Input variables:
%   AzVect: vector of azimuths
%   LinVect: vector of linear (X) variable 
%   testalfa: significance level for tests of significant correlation
%        tabled values are 0.10, 0.05, 0.01 for N < 100
%   fidO: output unit for external file (write if fidO>0)
% Output variable:
%   RtrnCode: return code

% Function called from this module:
%   Vector_Corr_CAssoc1Test.m
       
RtrnCode = 0;

% find nearest tabled value to specified alfa
alfalist = [0.10, 0.05, 0.01];
tblalfa = alfalist(1); diff = abs(testalfa - alfalist(1));
for iii = 2:3
    dd = abs(testalfa - alfalist(iii));
    if dd < diff
        diff = dd; tblalfa = alfalist(iii);
    end
end

NTot = length(AzVect);

     % get ranks in Theta and X
    
[xxx, Indx] = sort(AzVect);
jjj = 0;
for iii = 1:NTot
    jjj = jjj + 1;
    RankAz(Indx(iii)) = jjj;
end
%RankAz

[xxx, Indx] = sort(LinVect);
jjj = 0;
for iii = 1:NTot
    jjj = jjj + 1;
    RankLin(Indx(iii)) = jjj;
end
%RankLin

     % calculate measure of association and test significance
     
RankAz = 2*pi*RankAz/NTot;

Tc = RankLin*transpose(cos(RankAz));
Ts = RankLin*transpose(sin(RankAz));
Tc2Ts2 = Tc^2 + Ts^2;

if fix(NTot/2)*2 == NTot
    anxxx = cot(pi/NTot);
    an = 1/(1 + 5*anxxx^2 + 4*anxxx^4);
else
    an = (2*(sin(pi/NTot))^4)/(1 + cos(pi/NTot))^3;
end

Dn = an*Tc2Ts2;

Un = 24*Tc2Ts2/(NTot^2*(NTot + 1));
if NTot > 100
    alfa = testalfa;
else
    alfa = tblalfa;
end
cutoff = Vector_Corr_CAssoc1Test(NTot, alfa);
if NTot >= 40 
    Pval = chi2cdf(Un, 2);
end

     % output results
     
pt20 = sprintf('   Large-sample approximation to test distribution\n');       
pt21 = sprintf('     Estimated correlation (Dn) = %7.3g\n', Dn);
pt22 = sprintf('     Test statistic (Un) = %9.5g\n', Un);
pt22a = sprintf('     Significance level: alfa = %.3f\n', alfa);  
pt22B = sprintf('     No significance test: invalid alfa = %.5f\n',alfa);
pt23 = sprintf('     Test criterion (cutoff) = %9.5g\n', cutoff);
if Un > cutoff
    pt24 = sprintf('     Reject hypothesis of no association\n');
else
    pt24 = sprintf('     Cannot reject hypothesis of no association\n');
end
ptskip = sprintf('\n');

disp(pt20); disp(pt21); 
if cutoff > -990
    disp(pt22); disp(pt22a); disp(pt23); disp(pt24); 
else
    disp(pt22B)
end
if NTot >= 40
   pt25 = sprintf('     P-value = %7.4g\n', Pval);
   disp(pt25)
end
disp(ptskip)

if fidO > 0 
    fprintf(fidO, pt20); fprintf(fidO, pt21); 
    if cutoff > -990
        fprintf(fidO, pt22); fprintf(fidO, pt22a); fprintf(fidO, pt23);
        fprintf(fidO, pt24); 
    else
        fprintf(fidO, pt22B);
    end
    if NTot >= 40
       fprintf(fidO, pt25); 
    end
    fprintf(fidO, ptskip);
end

clear xxx Indx RankAz RankLin RArray DArray anxxx
