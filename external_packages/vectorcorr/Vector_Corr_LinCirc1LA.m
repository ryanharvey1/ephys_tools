function RtrnCode = Vector_Corr_LinCirc1LA(AzVect, LinVect, alfa, fidO)

% Vector_Corr_LinCirc1LA.m

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
% This uses the large-sample approximation to the test distribution
% Mardia and Jupp, 2000, p. 245 - 246
% Fisher, 1993, p. 145

% Input variables:
%   AzVect: vector of azimuths
%   LinVect: vector of linear (X) variable 
%   alfa: significance level for tests of significant correlation
%   fidO: output unit for external file (write if fidO>0)
% Output variable:
%   RtrnCode: return code
                
                                  
RtrnCode = 0;
                            
NTot = length(AzVect);

     % get single correlations to calculate multiple correlation
     
DArray = zeros(NTot, 3);
DArray(:,1) = cos(AzVect/57.3);
DArray(:,2) = sin(AzVect/57.3);
DArray(:,3) = LinVect;

RArray = corrcoef(DArray);

rcs = RArray(1,2);
rxs = RArray(1,3);
rxc = RArray(2,3);

R_AL2 = (rxc^2 + rxs^2 - 2*rxc*rxs*rcs)/(1 - rcs^2); 

     % test significance
     
U = ((NTot - 3)*R_AL2)/(1 - R_AL2);
cutoff = finv(1-alfa, 2, NTot-3);
Pval = 1 - fcdf(U, 2, NTot-3);

     % output results

ptskip = sprintf(' ');     
pt20A = sprintf('    Large-sample approximation to test distribution\n');
pt21 = sprintf('     Estimated correlation (R2xa) = %7.3g\n', R_AL2);
pt22 = sprintf('     Test statistic = %9.5g\n', U);
pt22a = sprintf('     Significance level: alfa = %.3f\n', alfa);  
pt23 = sprintf('     Test criterion (cutoff) = %9.5g\n', cutoff);
if U > cutoff
    pt24 = sprintf('     Reject hypothesis of no association\n');
else
    pt24 = sprintf('     Cannot reject hypothesis of no association\n');
end
pt25 = sprintf('     P-value = %6.4f\n', Pval);
pt26 = sprintf('     NOTE: Test and P-value approximate (small N)\n');

disp(pt20A); 
disp(pt21); disp(pt22); disp(pt22a); disp(pt23); disp(pt24); 
disp(pt25); 
if NTot < 30
    disp(pt26)
end
disp(ptskip)
if fidO > 0 
    fprintf(fidO, pt20A); 
    fprintf(fidO, pt21); fprintf(fidO, pt22); fprintf(fidO, pt22a);
    fprintf(fidO, pt23);
    fprintf(fidO, pt24); fprintf(fidO, pt25); 
    if NTot < 30
      fprintf(fidO, pt26);
    end
    fprintf(fidO, pt26);
end

clear xxx Indx RankAz RankLin RArray DArray anxxx