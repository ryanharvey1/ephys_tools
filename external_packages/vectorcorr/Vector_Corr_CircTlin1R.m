function RtrnCode = Vector_Corr_CircTlin1R(AzVect1, AzVect2, alfa, ...
                                          fidO, NTrials)

% Vector_Corr_CircTlin1R.m

% Copyright C 2005  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.
                                  
% circular-circular correlation  -  rhoT_hat
% T-linear association for constant shift between two vector variables
% Resampling used for test distribution
% Fisher, 1993, p. 151 - 153
% Mardia and Jupp, 2000, p. 250
% Resampling: Fisher, 1993, p. 214 - 218
% Also used for correlation of mean(theta) = f(linear) by first converting
%  linear to azimuth by AzL = 2 taninv(LV)
%  Fisher, p. 161

% Input variables:
%   AzVect1: vector of azimuths, variable 1
%   AzVect2: vector of azimuths, variable 2
%   alfa: significance level for tests of significant correlation
%   fidO: output unit for external file (write if fidO>0)
%   NTrials: number of resampling trials
% Output variable:
%   RtrnCode: return code

                                  
RtrnCode = 0;
                         
ptskip = sprintf('\n');
                                  
NTot = length(AzVect1);

% get basic calculations, including estimate of coefficient

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

% set up and loop over permutation trials

SaveStat = zeros(1, NTrials);
AzV3 = AzV2; AzV3 = 0;

for iii = 1:NTrials
        
    % get random permutation of linear variables
    
    R = randperm(NTot);
    
    for jjj = 1:NTot    
       AzV3(jjj) = AzV2(R(jjj));   
    end
     
     % calculate stat for measure-of-association coefficient

     cosAv2 = cos(AzV3);      sinAv2 = sin(AzV3);

     A = cosAv2*cosAv1;
     B = sinAv2*sinAv1;
     C = cosAv2*sinAv1;
     D = sinAv2*cosAv1;

     SaveStat(iii) = A*B - C*D;
     
 end
 
% sort calculated statistics and get cutoffs from generated distrib.

SortStat = sort(SaveStat);

L1 = fix(NTrials*alfa/2 + 0.5);
cutoff1 = SortStat(L1);
cutoff1 = 4*cutoff1 / sqrt((NTot^2 - E^2 - F^2)*(NTot^2 - G^2 - H^2));
L2 = NTrials - L1;
cutoff2 = SortStat(L2);
cutoff2 = 4*cutoff2 / sqrt((NTot^2 - E^2 - F^2)*(NTot^2 - G^2 - H^2));

     % output results
 
pt20 = sprintf(['   Test distribution based on resampling\n', ...
                '   Ref.: Fisher, 1993, p. 214 - 218\n']);        
pt21 = sprintf( '     Estimated correlation (rhoThat) = %7.3g\n', rhoT);
disp(pt20); disp(pt21)
if fidO > 0
    fprintf(fidO, pt20); fprintf(fidO, pt21); 
end
pt22alf = sprintf('     Significance level: alfa = %.3f\n', alfa);  
pt23 =    sprintf('     Test criteria (cutoffs): %7.3g,%7.3g\n', ...
                  cutoff1, cutoff2);
if rhoT < cutoff1 | rhoT > cutoff2 
   pt24 = sprintf('     Reject hypothesis of no association\n');
else
   pt24 = sprintf('     Cannot reject hypothesis of no association\n');
end
disp(pt22alf); disp(pt23); disp(pt24); 
disp(ptskip)
if fidO > 0 
      fprintf(fidO, pt22alf); fprintf(fidO, pt23); fprintf(fidO, pt24); 
      fprintf(fidO, ptskip);
end

clear x111 AzV1 AzV2 cosAv1 sinAv1 cosAv2 sinAv2 ...
      cos2Av1 sin2Av1 cos2Av2 sin2Av2 dAzV1 dAzV2 ...
      AzV3 iii jjj SaveStat SortStat R