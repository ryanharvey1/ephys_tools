function RtrnCode = Vector_Corr_PMCorr(LV, LV1, alfa, fidO)

% Vector_Corr_PMCorr.m

% Copyright C 2004  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.

% linear-linear  -  product-moment correlation, r

% Input variables:
%   LV: vector of data for variable 1
%   LV1: vector of data for variable 2
%   alfa: significance level for tests of significant correlation
%   fidO: output unit for external file (write if fidO>0)
% Output variable:
%   RtrnCode: return code

                                  
RtrnCode = 0;
                       
pt11 = sprintf('   Product-moment correlation\n');    
disp(pt11)
if fidO > 0 
   fprintf(fidO, pt11);
end

NTot = length(LV);

rcorr = corrcoef(LV, LV1);
r = zeros(1); r = rcorr(1,2);

     % test significance
     
tt = r*sqrt((NTot-2)/(1-r*r));
cutoff = tinv(1-alfa/2, NTot-2);

     % output results
     
pt21 = sprintf('     Estimated correlation (r) = %7.3g\n', r);
pt22 = sprintf('     Test statistic t = %8.3g\n', tt);
pt22a = sprintf('     Significance level: alfa = %.3f\n', alfa);  
pt23 = sprintf('     Test criterion (cutoff) = %7.3g\n', cutoff);
if abs(tt) > cutoff
    pt24 = sprintf('     Reject hypothesis of zero correlation\n');
else
    pt24 = sprintf('     Cannot reject hypothesis of zero correlation\n');
end
ptskip = sprintf('\n');

disp(pt21); disp(pt22); disp(pt22a); disp(pt23); disp(pt24); 
disp(ptskip)
if fidO > 0 
    fprintf(fidO, pt21); fprintf(fidO, pt22); fprintf(fidO, pt22a);
    fprintf(fidO, pt23); fprintf(fidO, pt24); 
    fprintf(fidO, ptskip);
end
