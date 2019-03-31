%DEMNPSI HOSA Toolbox Demo: Linear processes: Non-parametric system 
%			    identification;   Blind deconvolution. 
echo off

% demo of biceps, bicepsf and matul 

% A. Swami April 15, 1993, Nov 15 1994. 
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

% Linear processes: Non-parametric system identification / blind deconvolution

% Consider the linear system identification problem: the received
% signal is a linear process corrupted by additive noise, i.e., 
%
%   y(n) = sum {over i} h(i) u(n-i)   +  g(n) 
% 
% If we have access to input u(n), and some knowledge of the properties
% of noise g(n), we can use cross-correlation techniques to estimate
% h(n), the impulse response of the linear system. 

% However, in several problems such as seismic deconvolution, and 
% channel equalization (communications signals), we do not have 
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

% HOSA offers three routines for non-parametric (i.e., we do not assume
% an ARMA model) estimation of the impulse response, h(n). 
%
% BICEPS estimates the impulse response and the complex cepstrum, via 
%        third-order cumulants, and does not require phase unwrapping.
% BICEPSF is the frequency-domain version of BICEPS; for large orders
%	  (p and q, see below), this is faster. 
% MATUL  estimates the Fourier phase and magnitude of a signal using the
%        Matsuoka-Ulrych algorithm. 

load ma1 
% The test synthetic is a MA(3) synthetic with MA vector 
%  [1,0.9,0.385,-0.771];  the input sequence was i.i.d., and exponentially 
%  distributed;  white Gaussian noise was added to obtain a SNR of 20 dB.
% 
%  We will use biceps with p=8, q=8 

% Hit any key to continue
pause 

clf,
[he,ce] = biceps(y,8,8,128);
[val,loc] = max(abs(he));
n = length(he);   n = (n+1)/2;
set (gcf, 'Name', 'HOSA BICEPS')

disp(['Estimated IR has a max value of ',num2str(val), ...
      ' at i =',int2str(loc-n)])
% Recall that we always have a shift and scale ambiguity in this routine 
% Note that both the estimated cepstrum and impulse response die away rapidly 

% Hit any key to continue
pause 

% We will zoom in on the display 
iloc = abs(-n+loc)+3;
axis([-max(5,iloc), max(5,iloc), -val*1.1, val*1.1])
% We will shift and scale the true IR,  and display it by red circles 
htrue = [1, 0.9, 0.385, -0.771] * val;
hold on
plot((-n+loc):(-n+loc+3), htrue, 'ro'), 
hold off 
% Note that the estimated and true IR's are close to one another 
% In general, you will have to play with p and q parameters. 

% We can also use the FFT-based method BICEPSF, which is somewhat faster;
% there, one does not have to specify p or q, only the 
% the maximum number of third-order cumulants to be used. 

% Hit any key to continue
pause 
clc

% MATUL: Non-parametric magnitude and phase retrieval using the
%        Matsuoka-Ulrych algorithm. 
% 
%  We will compute the true third-order cumulants of an MA(2) model
%  with MA vector [1 -3.5 1.5]. 

% Hit any key to continue
pause

clf 
cmat = cumtrue([1 -3.5 1.5], 1, 3, 5); 
set (gcf, 'Name','HOSA MATUL') 

% We will compute the true bispectrum next: 

bsp = fft2(flipud(cmat),64,64);

% Note that we could have used routine BISPECT directly 
%
% Next, we use MATUL to estimate the impulse response 
hest = matul(bsp);

% 
% This routine should be used with care, since the phase unwrapping
% scheme is not very good.
% 
% Hit any key to return to the main menu 
pause 
echo off
clc

return


% Note that both the estimated cepstrum and impulse response die away rapidly 
% clf,
% [hf,cf] = bicepsf(y, 8);
% set (gcf, 'Name', 'HOSA BICEPSF')

