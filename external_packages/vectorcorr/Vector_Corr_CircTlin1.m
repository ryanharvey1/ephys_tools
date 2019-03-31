function RtrnCode = Vector_Corr_CircTlin1(AzVect1, AzVect2, alfa, ...
                                          fidO, flag, NTrials)

% Vector_Corr_CircTlin1.m

% Copyright C 2005  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.
                                  
% circular-circular correlation  -  rhoT_hat
% T-linear association for constant shift between two vector variables
% Uses large-sample approx. or resampling method
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
%   flag: 0=circular-circular   1=circular-linear converted to circ-circ
%   NTrials: number of resampling trials
% Output variable:
%   RtrnCode: return code

% Functions called from this module:
%   Vector_Corr_CircTlin1LA.m
%   Vector_Corr_CircTlin1R.m
                                
RtrnCode = 0;
                         
if flag == 0
   pt11 = sprintf(['   T-linear circular-circular association,\n', ...
                   '   assuming constant angular shift\n', ...
                   '   Refs.: Fisher, 1993, p. 151-153\n', ...   
                   '          Mardia and Jupp, 2000, p. 250\n']);   
else
   pt11 = sprintf(['   Linear-circular relation: ', ...
                   ' Mean(Theta) = function(X)\n',...
                   '   converted to circular-circular model\n', ...
                   '   Uses T-linear association, ', ...
                   'assuming constant angular shift\n', ...
                   '   Refs.: Fisher, 1993, p. 151-153, 161\n', ...   
                   '          Mardia and Jupp, 2000, p. 250\n']);   
end
           
disp(pt11)
if fidO > 0 
   fprintf(fidO, pt11);
end
ptlin = sprintf('     -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - \n');
ptskip = sprintf('\n');
                                  
NTot = length(AzVect1);

if NTot > 24
    
    % process assuming large-sample approximation -
    % overlapping resampling for N = 25 - 35
    
    RtrnCode = Vector_Corr_CircTlin1LA(AzVect1, AzVect2, alfa, fidO);
    
end

if NTot < 36
    
    % process assuming resampling (permutation test) -
    % overlapping large-sample approx. for N = 25 - 35
    
    RtrnCode = Vector_Corr_CircTlin1R(AzVect1, AzVect2, alfa, fidO,...
                                      NTrials);
    
end

disp(ptlin); disp(ptskip)
if fidO > 0 
      fprintf(fidO, ptlin); fprintf(fidO, ptskip);
end

