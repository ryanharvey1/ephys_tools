%DEMARMA  HOSA Toolbox Demo:  Linear Processes - Parametric (ARMA) model
%	  		      estimation;  Blind deconvolution.
echo off
% Demo of  arrcest, maest, armaqs, armarts,  aroder, maorder,
%	   cumtrue, bispect  (rpiid, armaysn, tls,   trench)


% A. Swami April 15, 1993; Nov 15, 1994
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

% Linear Processes - Parametric (ARMA) model estimation / Blind deconvolution
%
% Consider the linear system identification problem: the received
% signal is a linear process corrupted by additive noise, i.e.,
%
%   y(n) = sum {over i} h(i) u(n-i)   +  g(n)
%
% If we have access to input u(n), and some knowledge of the properties
% of noise g(n), we can use cross-correlation techniques to estimate
% h(n), the impulse response of the linear system.

% However, in several problems such as seismic deconvolution, and
% channel equalization of communications signals, we do not have
% access to the input u(n):  we must first estimate the impulse response
% h(n), so that we can estimate the input process.  This problem is
% known as blind deconvolution.  For the blind deconvolution problem
% to be well-posed, we must make some assumptions:  the usual assumption
% is that the input process is i.i.d..  If the process is also
% non-Gaussian, then, we can use cumulant-based techniques.  Recall
% that correlation based techniques cannot resolve the phase of the
% system (hence, the conventional assumption of minimum-phase), and
% that additive noise poses an additional problem.

% Hit any key to continue
pause
clc
% A parametric model is obtained by fitting an ARMA model to the
% the impulse response.
% An AR filter is an IIR filter with finite number of poles, and no zeros.
% A  MA filter is the same as a FIR filter (finite number of zeros, no poles).
% An ARMA filter is an IIR filter with finite number of zeros and poles.

% The HOSA Toolbox offers several routines for the estimation of the ARMA
%     parameters in the `blind deconvolution' problem.
%
% ARRCEST estimates AR parameters using correlation and cumulants
% MAEST   estimates MA parameters
% ARMAQS  estimates ARMA parameters via the q-slice algorithm
% ARMARTS estimates ARMA parameters via the residual time-series algorithm
% ARORDER estimates the order of an AR process, using correlation/cumulants
% MAORDER estimates the order of an MA process, using third-order cumulants
%
% Hit any key to continue
pause
%
%     Related routines:
% CUMTRUE computes the theoretical cumulants of an ARMA model
% BISPECT computes the theoretical bispectrum of an ARMA model
% RPIID   generates sequences of i.i.d. random variables
% ARMASYN generates ARMA synthetics
% TLS     solves a linear system of equations via Total Least Squares
% TRENCH  Trench recursion for non-symmetric Toeplitz matrices

% Hit any key to continue
pause
clc
echo off

l_arma = str2mat( ...
        'AR parameter estimation    ',  ...
	'MA parameter estimation    ',  ...
	'ARMA parameter estimation  ', ...
	'Model order determination  ', ...
	'True cumulants/polyspectra ');

c_arma = str2mat('d_ar','d_ma','d_arma','d_order','d_true');
choices ('HosatArmaDemo','HOSA - Parametric estimators', ...
	  l_arma, c_arma, 1);

echo off
clc
return