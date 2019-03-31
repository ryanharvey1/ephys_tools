%DEMADAPT HOSA Toolbox Demo: Adaptive AR parameter estimation
%
echo off
% demo of RIVDL and RIVTR (which call IVCAL)

% A. Swami April 15, 1993,  Nov 15, 1994.
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

%		 Adaptive AR parameter estimation

% The problem of obtaining adaptive estimates of the parameters of an AR
% process occurs in several situations: for example, in communication signals
% (adaptive channel equalization), in seismic (adaptive deconvolution, to
% handle the non-stationarity of the signal), and, of course, in situations
% where the AR parameters themselves change with time (e.g., estimation of
% the instantaneous frequency of a chirp signal).

% The HOSA Toolbox offers two routines for adaptive linear prediction:
% 1. RIVDL : double lattice form of the recursive instrumental variable method
% 2. RIVTR : transversal form of the recursive instrumental variable method
%
%        (IVCAL is used to compute the instrumental variables.)

%Hit any key to continue
pause

% We will test the algorithms on noise-free and noisy data
% A zero-mean, i.i.d., exponentially distributed sequence was convolved
%   with the AR filter [1,-1.5,0.8], to generate the signal "y".
% The signal was corrupted by white Gaussian noise to obtain  "zw",
%   which has a SNR of 10 dB.
% The signal was corrupted by colored Gaussian noise to obtain "zc",
%   which has a SNR of 10 dB;  colored noise was obtained by passing
%   white noise through the AR filter [1,0,0.49] .

%Hit any key to continue
pause
echo off

l_adapt = str2mat('Lattice filters - rivdl', ...
	         'Transversal filters - rivtr');
c_adapt = str2mat('d_rivdl','d_rivtr');

choices('HosatAdaptDemo','HOSA - Adaptive estimators', l_adapt, c_adapt,1);
echo off
clc
return