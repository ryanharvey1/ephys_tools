function RtrnCode = Vector_Corr_CircTlin1LA(AzVect1, AzVect2, alfa, ...
                                          fidO)

% Vector_Corr_CircTlin1LA.m

% Copyright C 2005  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.
                                  
% circular-circular correlation  -  rhoT_hat
% T-linear association for constant shift between two vector variables
% Fisher, 1993, p. 151 - 153
% Mardia and Jupp, 2000, p. 250
% Uses large-sample approximation
% Also used for correlation of mean(theta) = f(linear) by first converting
%  linear to azimuth by AzL = 2 taninv(LV)
%  Fisher, p. 161

% Input variables:
%   AzVect1: vector of azimuths, variable 1
%   AzVect2: vector of azimuths, variable 2
%   alfa: significance level for tests of significant correlation
%   fidO: output unit for external file (write if fidO>0)
% Output variable:
%   RtrnCode: return code

% Functions called from this module:
%   VectMean_arctan.m

                                  
RtrnCode = 0;              
                                  
NTot = length(AzVect1);

     % compute sums to be used in calculations
     
AzV1 = AzVect1/57.3;     AzV2 = AzVect2/57.3;
cosAv1 = cos(AzV1);      sinAv1 = sin(AzV1);
cosAv2 = cos(AzV2);      sinAv2 = sin(AzV2);
cos2Av1 = cos(2*AzV1);   sin2Av1 = sin(2*AzV1);
cos2Av2 = cos(2*AzV2);   sin2Av2 = sin(2*AzV2);

A = transpose(cosAv1)*cosAv2;
B = transpose(sinAv1)*sinAv2;
C = transpose(cosAv1)*sinAv2;
D = transpose(sinAv1)*cosAv2;
x111 = ones(1,NTot);
E = x111*cos2Av1;
F = x111*sin2Av1;
G = x111*cos2Av2;
H = x111*sin2Av2;

     % compute rhoT as estimate of correlation

rhoT = 4*(A*B - C*D);
rhoT = rhoT / sqrt((NTot^2 - E^2 - F^2)*(NTot^2 - G^2 - H^2));

     % output estimate

ptskip = sprintf('\n');
pt20 = sprintf('   Large-sample approximation to test distribution\n');     
pt21 = sprintf('     Estimated correlation (rhoThat) = %7.3g\n', rhoT);
disp(pt20); disp(pt21)
if fidO > 0
    fprintf(fidO, pt20); fprintf(fidO, pt21); 
end

     % set up for making tests of significance
     % Fisher, 1993, p. 152

    UU1 = x111*cosAv1;      VV1 = x111*sinAv1;
    UU2 = x111*cosAv2;      VV2 = x111*sinAv2;
    R1 = sqrt(UU1^2 + VV1^2)/NTot;
    R2 = sqrt(UU2^2 + VV2^2)/NTot;
    Thet1 = VectMean_arctan(VV1, UU1)/57.3;
    Thet2 = VectMean_arctan(VV2, UU2)/57.3;
    
    dAzV1 = 2*(AzV1 - Thet1);       dAzV2 = 2*(AzV2 - Thet2);
    cos2Av1 = cos(dAzV1);           sin2Av1 = sin(dAzV1);
    cos2Av2 = cos(dAzV2);           sin2Av2 = sin(dAzV2);
    alf1 = x111*cos2Av1/NTot;       alf2 = x111*cos2Av2/NTot;  
    bet1 = x111*sin2Av1/NTot;       bet2 = x111*sin2Av2/NTot;
    
    U1 = (1 - alf1^2 - bet1^2)/2;   U2 = (1 - alf2^2 - bet2^2)/2;
    V1 = R1^2*(1 - alf1);           V2 = R2^2*(1 - alf2);

     % test significance

   % case i
   pt51 = sprintf(['     Assume either distribution has vector-length',...
                   ' R = 0\n']);
   U = abs(NTot*rhoT);         
   cutoff = -log(alfa);
   pt22 = sprintf('         Test statistic = %9.5g\n', U);
   pt22a = sprintf('         Significance level: alfa = %.3f\n', alfa);  
   pt23 = sprintf('         Test criterion (cutoff) = %9.5g\n', cutoff);
   if U > cutoff
     pt24 = sprintf('         Reject hypothesis of no association\n');
   else
     pt24=sprintf('         Cannot reject hypothesis of no association\n');
   end
   Pval = exp(-U);
   pt25 = sprintf('         P-value = %7.4g\n', Pval);
   disp(pt51); disp(pt22); disp(pt22a); disp(pt23); disp(pt24); disp(pt25)
   if NTot < 35
       pt26 = sprintf('      NOTE: Test and P-value approximate (small N)\n');
       disp(pt26)
   end
   if fidO > 0 
     fprintf(fidO, pt51); fprintf(fidO, pt22); fprintf(fidO, pt22a);
     fprintf(fidO, pt23);
     fprintf(fidO, pt24); fprintf(fidO, pt25);
     if NTot < 35
       disp(pt26)
     end
   end

   % case ii
   pt51 = sprintf(['     No assumption of zero-length vector in ', ...
                   'either distribution\n']);
   ZZZ = sqrt(NTot)*U1*U2*rhoT/sqrt(V1*V2);
   cutoff = norminv(1 - alfa/2);
   pt22 = sprintf('         Test statistic (Z) = %9.5g\n', ZZZ);
   pt23 = sprintf('         Test criterion (cutoff) = %9.5g\n', cutoff);
   if abs(ZZZ) > cutoff
     pt24 = sprintf('         Reject hypothesis of no association\n');
   else
     pt24=sprintf('         Cannot reject hypothesis of no association\n');
   end
   Pval = 2*(1 - normcdf(abs(ZZZ),0,1));
   pt25 = sprintf('         P-value = %7.4g\n', Pval);
   disp(pt51); disp(pt22); disp(pt22a); disp(pt23); disp(pt24); disp(pt25)
   if NTot < 35
       disp(pt26)
   end
   disp(ptskip)
   if fidO > 0 
      fprintf(fidO, pt51); fprintf(fidO, pt22); fprintf(fidO, pt22a);
      fprintf(fidO, pt23);
      fprintf(fidO, pt24); fprintf(fidO, pt25);
      if NTot < 35
          fprintf(fidO, pt26)
      end
      fprintf(fidO, ptskip);
   end

clear x111 AzV1 AzV2 cosAv1 sinAv1 cosAv2 sinAv2 ...
      cos2Av1 sin2Av1 cos2Av2 sin2Av2 dAzV1 dAzV2