%D_BISP  HOSA Demo: Estimating the Bispectrum (Direct & Indirect methods)
%	             and Bicoherence (auto and cross)
% demos of  bispecd, bispeci, bicoher, bispecdx, bicoherx
%       
echo off
% A. Swami Oct 18, 1997.
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
%        Estimating Cumulants, Bispectra and Bicoherence

%  CUMEST   - estimates second-, third-, and fourth-order cumulants.
%  CUM2X, CUM3X, CUM4X - estimate 2nd, 3rd and 4th order cross-cumulants
%  BISPECD  - direct (FFT-based) estimate of the  auto-bispectrum
%  BISPECI  - indirect estimate of the auto-bispectrum
%  BICOHER  - direct (FFT-based) estimate of the bicoherence
%  BISPECDX - direct (FFT-based) estimate of the cross-bispectrum
%  BICOHERX - direct (FFT-based) estimate of the cross-bicoherence
%
% The power spectrum is defined as the Fourier Transform of the
% autocorrelation.   Similarly, the bispectrum is defined as the 2-D Fourier
% transform of the third-order cumulant C3(i,j). The Fourier transform of the
% sample estimate of the auto-correlation does not yield a consistent (i.e.,
% reliable) estimate of the power spectrum - the estimate must be smoothed,
% either by using a lag window or smoothing in the frequency domain.
% Similarly, the Fourier transform of the third-order cumulant, by itself,
% does not yield a consistent (reliable) estimate of the bispectrum.  The
% estimate must be smoothed, either in the time domain or in the frequency
% domain.   This leads to two different methods of estimating the bispectrum.

% Hit any key to continue
pause
clc

echo off
l_bisp = str2mat( ...
	 'Indirect estimate (bispeci) ', ...
         'Direct estimate   (bispecd) ', ...
	 'Cross-bispectrum  (bispecdx) ', ...
	 'Auto-bicoherence  (bicoher)  ', ...
         'Cross-bicoherence (bicoherx) ' );

c_bisp = str2mat( ...
	   'd_bispi', 'd_bispd', 'd_bispdx', 'd_bic', 'd_bicx' );

choices ('HosatBispDemo', 'HOSA - bispectra and bicoherence', ...
	   l_bisp, c_bisp, 1);

echo off
return