%D_TRUE  HOSA Demo for computing true cumulants / polyspectra
%	demo of cumtrue, bispect, trispect

echo off

% A. Swami Jan 20, 1995
% Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.6 $

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

%	   Computing true cumulants and polyspectra
%
% HOSA offfers the following routines
%
% CUMTRUE  - computes true cumulants of ARMA models
% BISPECT  - computes the theoretical bispectrum of an ARMA process
% TRISPECT - computes 1-D slices of the theoretical trispectrum of an
%		ARMA process
%
% Related Routines:
% TRENCH: Estimates the AR parameters and the reflection coefficients given
%         a slice of a cumulant sequence, where the corresponding matrix
%         is Toeplitz, but non-symmetric.

%Hit any key to continue
pause

echo off

l_true = str2mat('True cumulants', 'True bispectrum', 'True trispectrum');
c_true = str2mat('d_truec', 'd_trueb', 'd_truet');
choices('HosatTrueDemo','Theoretical polyspectra',l_true,c_true,1);

%

% Hit any key to return to main menu
pause
echo off
clc
return