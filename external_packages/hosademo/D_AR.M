%D_AR   HOSA Demo:  Linear Processes - Parametric (AR) model estimation
%      Demo of arrcest

echo off
% A. Swami Oct 18, 1997.
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

%              AR parameter estimation methods: 
%
% These routines use second- and/or third- and/or fourth-order cumulants
% to fit AR models to non-Gaussian processes.  Routine ARRCEST uses the
% so-called 'normal' equations based on cumulants and/or correlation,
% to estimate the AR parameters. 
% 
% The test synthetic "y" is an AR(2) synthetic;  the true AR parameters
% were [1 -1.5 0.8]; input was i.i.d. and exponentially distributed. 
% Additive white Gaussian noise was added to obtain a SNR of 20 dB. 
%
% We will use ARRCEST to estimate the AR parameters

% Hit any key to continue
pause 
load ar1 
ar(:,1) = arrcest(y,2,0, 2,12,128); 
ar(:,2) = arrcest(y,2,0, 3,12,128); 
ar(:,3) = arrcest(y,2,0, 4,12,128); 
ar(:,4) = arrcest(y,2,0,-3,12,128); 
ar(:,5) = arrcest(y,2,0,-4,12,128); 

%       1         2         3         4         5 
% ---------------------------------------------------------- 
disp(ar) 
% 
% The five columns correspond to five different estimates using
% 1. correlation 
% 2. third-order cumulants
% 3. fourth-order cumulants
% 4. correlation and third-order cumulants
% 5. correlation and fourth-order cumulants
%
% Since the process is linear non-Gaussian, at high SNR, the various
% estimates are close to one another. 

% Hit any key to continue
pause 
echo off 
clc

