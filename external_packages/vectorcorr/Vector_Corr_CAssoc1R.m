function RtrnCode = Vector_Corr_CAssoc1R(AzVect, LinVect, testalfa, ...
                                        fidO, NTrials)

% Vector_Corr_CAssoc1R.m

% Copyright C 2005  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.

% Linear-vector association   -   Un, Dn                               
% Distribution-free rank method - watch out: sensitive to tied data
% Uses resampling (permutation test) for distribution of test
% Mardia and Jupp, 2000, p. 246-248
% Fisher, 1993, p. 140-141
% Resampling: Fisher, 1993, p. 214 - 218

% Input variables:
%   AzVect: vector of azimuths
%   LinVect: vector of linear (X) variable 
%   testalfa: significance level for tests of significant correlation
%   fidO: output unit for external file (write if fidO>0)
%   NTrials: number of resampling trials
% Output variable:
%   RtrnCode: return code

       
RtrnCode = 0;

NTot = length(AzVect);

% get basic calculations, including estimate of coefficient

     % get ranks in Theta and X
    
[xxx, Indx] = sort(AzVect);
jjj = 0;
for iii = 1:NTot
    jjj = jjj + 1;
    RankAz(Indx(iii)) = jjj;
end

RankAz = 2*pi*RankAz/NTot;

[xxx, Indx] = sort(LinVect);
jjj = 0;
for iii = 1:NTot
    jjj = jjj + 1;
    RankLin(Indx(iii)) = jjj;
end

     % calculate measure of association 

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

% set up and loop over permutation trials

SaveStat = zeros(1, NTrials);
LinVect2 = zeros(1, NTot);

for iii = 1:NTrials
        
    % get random permutation of linear variables
    
    R = randperm(NTot);
    
    for jjj = 1:NTot    
       LinVect2(jjj) = LinVect(R(jjj));   
    end
     
     % calculate stat for measure-of-association coefficient

    [xxx, Indx] = sort(LinVect2);
    ijjj = 0;
    for iiii = 1:NTot
        ijjj = ijjj + 1;
        RankLin(Indx(iiii)) = ijjj;
    end

    Tc = RankLin*transpose(cos(RankAz));
    Ts = RankLin*transpose(sin(RankAz));
    Tc2Ts2 = Tc^2 + Ts^2;

    SaveStat(iii) = Tc2Ts2;
    
end

% sort calculated statistics and get cutoffs from generated distrib.

SortStat = sort(SaveStat);

L1 = fix(NTrials*(1 - testalfa) + 0.5);
cutoff = an*SortStat(L1);

alfa = testalfa;

     % output results

pt20 = sprintf(['   Test distribution is based on resampling\n', ...
               '   Ref.: Fisher, 1993, p. 214 - 218\n']);  
pt21 = sprintf('     Estimated correlation (Dn) = %7.3g\n', Dn);
pt22a = sprintf('     Significance level: alfa = %.3f\n', alfa);  
pt23 = sprintf('     Test criterion (cutoff) = %7.3g\n', cutoff);
if Dn > cutoff
    pt24 = sprintf('     Reject hypothesis of no association\n');
else
    pt24 = sprintf('     Cannot reject hypothesis of no association\n');
end

ptskip = sprintf('\n');

disp(pt20); disp(pt21); 
disp(pt22a); disp(pt23); disp(pt24); 
disp(ptskip)
if fidO > 0 
    fprintf(fidO, pt20); fprintf(fidO, pt21); 
    fprintf(fidO, pt22a); fprintf(fidO, pt23);
    fprintf(fidO, pt24); 
    fprintf(fidO, ptskip);
end

clear xxx Indx RankAz RankLin anxxx SaveStat SortStat R iii jjj LinVect2
clear ijjj iiii 
