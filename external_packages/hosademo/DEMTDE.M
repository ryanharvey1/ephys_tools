%D_TDELAY  HOSA Demo: Time Delay Estimation (TDE)
echo off

% demos of tde, tdeb, tder  (tdegen)

% A. Swami April 15, 1993
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

clear , clc,
echo on

%                 Time Delay Estimation (TDE)
%
%  The TDE problem occurs in several real-world situations, such as
%  sonar     (estimate the range and bearing of an underwater acoustic source)
%  radar       (estimate the range and bearing of a radio transmitter)
%  geophysics  (estimate the epicenter of an earthquake)
%  biomedicine (estimate the location of `dipoles' in the brain)
%  material sciences (estimate some intrinsic property of an object,
%                     by measuring the flight time of a signal)
%
%  The signal is recorded at two or more sensors located in a known spatial
%  configuration;  estimates of the delay in arrival time of the signal
%  from sensor to sensor are sought;  these time delay estimates
%  are then used to estimate the source and bearing (for the source
%  localization problem) or to infer some properties of the medium of
%  propagation.

% Hit any key to continue
pause
clc
%  In practice, the sensor measurements are contaminated by noise which
%  may be spatially correlated.  If the spatial correlation matrix of
%  the sensor noises is known, then, second-order techniques may be used
%  (pre-whiten and cross-correlate) to estimate the signal delay.  If
%  the noise correlation matrix is not known, we cannot correctly pre-whiten
%  the data: if we cross-correlate the signals, we may not be able to
%  differentiate between the signal delay and the noise delay:
%  second-order statistics are inadequate.

%  If the signals are non-Gaussian, or deterministic, and the noises
%  are Gaussian, we can use third- or higher-order cumulants to estimate
%  the time delays.  If the signal has non-zero third order cumulants,
%  then, the noises may be symmetric distributed (e.g., Gaussian, Laplace,
%  Uniform, etc.).

% Hit any key to continue
pause
clc

% HOSA offers the following routines related to TDE
%  TDEGEN - generates synthetics for the two-sensor TDE problem
%  TDE    - estimates the time delay using cumulants
%  TDER   - estimates the time delay using the ML-windowed cross-correlation
%  TDEB   - estimates the time delay using the bispectrum
%
%   Related routines are DOA and DOAGEN (bearing estimation, direction of
%   arrival)
%
% Hit any key to continue
pause
clc
echo off
l_tde = str2mat( ...
	'Cross-cumulant method (tde) ', ...
        'Cross-bispectrum method (tdeb) ', ...
	'ML-window cross-correlation method (tder) ' ) ;
c_tde = str2mat('d_tde','d_tdeb','d_tder');
choices ('HosatTdeDemo',' HOSA - Time Delay Estimation',l_tde, c_tde, 1);
echo off
clc
return

%
% 1. TDEGEN may be used to generate synthetics for the TDE problem.
%    A signal and its delayed version, contaminated by noise sequences
%    which are correlated with each other, are generated.
%


%
% In this example, the parametric method (TDE) was better than the
% hologram method (TDEB), which in turn was better than the ML-correlation
% method (TDER).
%