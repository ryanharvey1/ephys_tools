function  testValue = Vector_Corr_CirRankDeltaTest(N, alfa, flag)

%Vector_Corr_CirRankDeltaTest.m

% Copyright C 2004  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.

% Look up cutoff values for DELTAhatn test for testing 
% circular-circular rank correlation for independence.
% Tables referenced below.

% NOTE: Not executed by Vector_Corr
%       Reference on method to calculate DELTAhatn
%       Fisher, 1993, p. 146 - 148

% Variables input:
%   N: sample size
%   alfa: signficance level to use
% Variable output:
%   testValue:  Cutoff value for test

% -----------------------------------------------------------
  
    % test for DELTA hat n

%TableR contains N in col.1, cutoffs for various alfa in col.2-7
%  Alfa values are 0.10, 0.05, 0.025, 0.01, 0.005, 0.001 
%  (from Fisher, 1993, p. 235)

TableR = [8 2.10 2.78 3.40 4.25 4.74 6.23;
          9 2.07 2.72 3.33 4.16 4.64 6.08;
         10 2.04 2.68 3.28 4.09 4.56 5.96;
         11 2.01 2.65 3.24 4.02 4.50 5.85;
         12 1.99 2.62 3.20 3.97 4.44 5.77;
         13 1.97 2.59 3.17 3.93 4.40 5.70;
         14 1.96 2.57 3.14 3.89 4.36 5.64;
         15 1.95 2.56 3.12 3.86 4.33 5.59;
         20 1.90 2.49 3.04 3.75 4.21 5.40;
         25 1.87 2.46 2.99 3.68 4.14 5.29;
         30 1.86 2.43 2.96 3.64 4.09 5.22;
        200 1.77 2.31 2.81 3.42 3.85 4.85; 
       9999 1.77 2.31 2.81 3.42 3.85 4.85]; 
    
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
           else if alfa == 0.001
                   alfacol = 7;
               end
           end
           end
       end
   end
   end
   
   if alfacol == 0
       return
   end

   % get and return cutoff value for given N and alfa        

   testValue = interp1(TableR(:,1), TableR(:,alfacol), N, 'cubic'); 

