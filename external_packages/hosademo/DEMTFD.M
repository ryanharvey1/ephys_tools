%D_TFD   HOSA Demo of Time-Frequency Distributions:
%        second, third and fourth-order Wigner distributions.
%
echo off
% demo of Wigner routines (wig2,wig2c, wig3,wig3c, wig4,wig4c)

% A. Swami Nov 15, 1994
% Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.8 $

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

%		   Time-Frequency Distributions:
%        second, third and fourth-order Wigner distributions.

%  The power spectrum, the bispectrum and the trispectrum are appropriate
%  tools to study stationary signals.  In the case of non-statioanry signals,
%  or overlapping deterministic signals,  time-frequency distributions (TFD's)
%  such as those provided by the conventional and higher-order Wigner
%  transforms are useful.  These TFD's essentially compute the Fourier
%  Transform of an appropriately defined instantaneous estimate of the
%  moment function.
%
%  HOSA offers several routines:
%  WIG2 and WIG2C - conventional or second-order Wigner spectrum
%  WIG3 and WIG3C - diagonal slice of the Wigner bispectrum
%  WIG4 and WIG4C - diagonal slice of the Wigner trispectrum
%   The suffix C indicates that a generalized form of the Choi-Williams filter
%   is used in order to suppress the effect of cross terms; it should be noted
%   that such a filter used in conjunction with the Wigner bispectrum
%   may destroy the signal; hence, routine WIGC should be used with care.
%
% Hit any key to continue
pause
clc

% First let us look at the test signals
load wigdat
clf
subplot(211), plot(s2),  title('s2'), grid on
subplot(212), plot(s3),  title('s3'), grid on
set(gcf, 'Name','HOSA WIGDAT')

% Signal s2 is a harmonic with a nominal frequency of 0.05 Hz, which has
% been multiplied by a Gaussian window, and is centered at n=20.
%
% Signal s3 is a harmonic with a nominal frequency of 0.15 Hz, which has
% been multiplied by a Gaussian window, and is centered at n=50.

% Hit any key to continue
pause
clc

echo off
l_wig = str2mat('Wigner Spectrum ', ...
		'Wigner Bispectrum ', ...
		'Wigner Trispectrum ');

c_wig = str2mat('d_wig2', 'd_wig3', 'd_wig4' );
choices('HosatTfdDemo','HOSA - Time-frequency distributions', l_wig,c_wig,1);

echo off
clc
return