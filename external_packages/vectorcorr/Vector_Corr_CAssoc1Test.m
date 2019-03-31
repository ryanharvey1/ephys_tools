function  UnValue = Vector_Corr_CAssoc1Test(N, alfa)

%Vector_Cor_CAssoc1Test.m

% Copyright C 2004  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.

%Look up cutoff values for C Association test (Un) for
%testing X vs theta independence  

% Variables input:
%   N: sample size
%   alfa: significance level to use 
%         Tabled values (0.10, 0.05, 0.01)
% Variable output:
%   UnValue:  Cutoff value for Un test of C Association
%             -999 implies no test due to incorrect alfa

% -----------------------------------------------------------

%TableUn contains N in col.1, cutoffs for various alfa in col.2-4
%  Alfa values are 0.10, 0.05, 0.01 
%  (from Mardia and Jupp, 2000, p. 380)

TableUn= [6  4.57  4.67  4.95;
          7  4.30  4.90  5.75;
          8  4.49  5.17  6.15;
          9  4.50  5.34  6.68;
         10  4.52  5.48  6.68;
         11  4.55  5.5   7.2;
         12  4.57  5.6   7.5;
         15  4.59  5.7   7.9;
         20  4.60  5.8   8.3;
         30  4.60  5.9   8.7;
         40  4.60  5.9   8.8;
         50  4.61  6.0   8.9;
        100  4.61  6.0   9.1]; 
    
   UnValue = -999;

   % determine which alfa column to use
   
   alfacol = 0;
   if alfa == 0.10
       alfacol = 2;
   else if alfa == 0.05
       alfacol = 3;
   else if alfa == 0.01 
       alfacol = 4;
       end
   end
   end

   % get and return cutoff value for given N and alfa        

   if N <= 100   
       if alfacol == 0 
          return
       end   
       UnValue = interp1(TableUn(:,1), TableUn(:,alfacol), N, 'cubic'); 
   else
       UnValue = chi2inv(1-alfa, 2); 
   end  

     