%D_MAORD  HOSA Demo:  Linear Processes - MA model order determination.
echo off
% Demo of  maorder

% A. Swami Jan 20, 1995 
% Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.5 $

%     RESTRICTED RIGHTS LEGEND
% Use, duplication, or disclosure by the Government is subject to
% restrictions as set forth in subparagraph (c) (1) (ii) of the 
% Rights in Technical Data and Computer Software clause of DFARS
% 252.227-7013. 
% Manufacturer: United Signals & Systems, Inc., P.O. Box 2374, 
% Culver City, California 90231. 
%
%  This material may be reproduced by or for the U.S. Government pursuant 
%  to the copyright license under the clause at DFARS 252.227-7013. 

clear, clc, 

echo on

%                 MAORDER: MA order determination
%
% We will use MAORDER to determine the order of an MA(3) process:
% the true MA parameters are [1, 0.9, 0.385, -0.771];  input was i.i.d. 
% exponential; white Gaussian noise was added to obtain an SNR of 20 dB
% We will use the default PFA (probability of false alarm) value of 0.05 

%Hit any key to continue
pause  
 load ma1
 q = maorder (y, 0, 6); 

% If the absolute value of the estimated C(q,0) exceeds the threshold, then
% the last column will be 0, indicating that the data are not consistent with
% a MA(q) model.  From the table we conclude that the data are inconsistent
% with  MA(0) and MA(2) models; they are consistent with model orders of
% 3, 4, 5 and 6; hence, we choose the order 3; of course, the true order may
% be larger than the specified maximum order of 6. 

% Hit any key to continue
pause
echo off 
clc 

