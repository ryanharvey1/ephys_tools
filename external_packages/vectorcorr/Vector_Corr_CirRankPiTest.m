function  testValue = Vector_Corr_CirRankPiTest(N, alfa)

%Vector_Corr_CirRankPiTest.m

% Copyright C 2004  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.

% Look up cutoff values for PIhatn test for testing 
% circular-circular rank correlation for independence.
% Tables referenced below.

% Variables input:
%   N: sample size
%   alfa: significance level to use
% Variable output:
%   testValue:  Cutoff value for test

% -----------------------------------------------------------
    
    % PIhat n test
 
%TableC contains N in col.1, cutoffs for various alfa in col.2-6
%  Alfa values are 0.10, 0.05, 0.025, 0.01, 0.005
%  (from Fisher, 1993, p. 237)

TableC = [8 1.80 2.55 3.26 4.11 4.80;
          9 1.78 2.52 3.23 4.09 4.78;
         10 1.76 2.50 3.21 4.07 4.77;
         11 1.74 2.48 3.19 4.06 4.75;
         12 1.73 2.46 3.17 4.04 4.74;
         13 1.72 2.45 3.15 4.03 4.73;
         14 1.71 2.44 3.14 4.02 4.72;
         15 1.71 2.43 3.13 4.02 4.71;
         20 1.68 2.40 3.10 3.99 4.67;
         25 1.67 2.38 3.07 3.97 4.64;
         30 1.66 2.36 3.06 3.96 4.63;
        200 1.61 2.30 2.99 3.91 4.60; 
       9999 1.61 2.30 2.99 3.91 4.60]; 
     
   testValue = -999;

   % determine which alfa column to use
   
   alfacol = 0;
   if alfa == 0.10
       alfacol = 2;
   else if alfa == 0.05
       alfacol = 3;
   else if alfa == 0.025
           alfacol = 4;
       else if alfa == 0.01
           alfacol = 5;
       else if alfa == 0.005
               alfacol = 6;
           end
       end
       end
   end
   end
   
   if alfacol == 0
       return
   end

   % get and return cutoff value for given N and alfa        

   testValue = interp1(TableC(:,1), TableC(:,alfacol), N, 'cubic'); 
 
     