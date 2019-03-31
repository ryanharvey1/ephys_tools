%D_ARMA  HOSA Demo:  Linear Processes - Parametric (ARMA) model estimation
%                        Blind deconvolution. 
echo off
% Demo of  armaqs, armarts

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


%             ARMA parameter estimation methods 
%
% HOSA offers routines, ARMAQS and ARMARTS to estimate
% the parameters of ARMA models. 
%
% The test synthetic "z" is an ARMA(2,1) synthetic with AR parameters
% [1,-0.8,0.65] and MA parameters [1,-2]; input was i.i.d. and exponentially 
% distributed; white Gaussian noise was added to obtain a SNR of 20 dB. 
% 
% We will use ARMAQS to estimate the parameters: 
% This implements the `q-slice' method which obtains simultaneous 
% estimates of the AR and MA parameters. 
% We will use p=2, q=1 (true model orders), and 3-rd order cumulants. 

%Hit any key to continue
pause

load arma1

[avec,bvec] = armaqs(y,2,1,3,10,128);
%  AR estimate: 
disp(avec'),
%  MA estimate: 
disp(bvec')

% Hit any key to continue
pause

% Next, we will use ARMARTS (residual time-series method)
% This algorithm uses a three step procedure; it uses ARRCEST to esimate the
% AR parameters; next, it deconvolves the AR part, and then uses MAEST to
% estimate the MA parameters. 
% We will use p=2, q=1 (true model orders), and 3-rd order cumulants. 

%Hit any key to continue
pause


[avec,bvec] = armarts(y,2,1,3,12,128);
%  AR estimate: 
disp(avec'),
%  MA estimate: 
disp(bvec')

% Hit any key to continue
pause 
echo off
clc
