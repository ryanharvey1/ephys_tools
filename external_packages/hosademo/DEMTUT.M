%DEMTUT	  A quick tutorial introduction to Higher-Order Statistics 

echo off 

% A. Swami April 15, 1993, Nov 15, 1994. 
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

%	 What are higher-order statistics ? 
% By higher-order statistics, we mean the higher-order moments and
% certain non-linear combinations of the moments, called CUMULANTS.

% Let x(n) denote a discrete-time, real-valued stationary random process.
% We will assume that all the statistical quantities defined below,
% exist, and are finite-valued. 
% The second-order moment function is defined by 
%
%   M2(i) = E { x(n) x(n+i) } 
%
% where E denotes the statistical operation, and M2 is a mnemonic 
% for second-order moment.
% This is the familiar autocorrelation function. 
 
% Hit any key to continue 
pause
clc

% The higher-order moment functions are defined as the expected values
% of multiple products: they are multi-dimensional sequences.
% For example, the third-order moment is the expected value of a 
% triple product, 
%
%    M3(i,j) = E { x(n) x(n+i) x(n+j) } 
%
% Note that M3 is a two-dimensional sequence. 
% Similarly, the fourth-order moment is the expected value of the 
% product of four terms, namely, 
%
%    M4(i,j,k) = E { x(n) x(n+i) x(n+j) x(n+k) } 
%
% Note that M4 is a three-dimensional sequence. 

% Hit any key to continue
pause 
clc

% Cumulants are non-linear combinations of the moments of the process.
% Unlike moments, the cumulants of a process are invariant to changes 
% in the mean value of the process. 
% Hence, let us make the convenient assumption that the process has zero mean.
% The second-order cumulant, denoted by C2(i), is the auto-covariance 
% function.  
% For a zero-mean process, the third-order cumulant is also identical
% with the third-order moment of the process, i.e., 
% 
%  C2(i)   = M2(i)   = E { x(n) x(n+i) }
%  C3(i,j) = M3(i,j) = E { x(n) x(n+j) x(n+k) } 

% Things start getting a bit messy for orders larger than three
% Thus, the fourth-order cumulant is defined via 
%
% C4(i,j,k) = M4(i,j,k) - M2(i) * M2(j-k) - M2(j) * M2(k-i) - M2(k) * M2(i-j) 
%
%  Why should we worry about these multi-dimensional statistics, 
%  and their non-linear combinations ? 

% Hit any key to continue 
pause 
clc

% If x(n) is a Gaussian process, then we know that its statistics are
% completely characterized by its autocorrelation function: 
% all of its higher-order moments can be defined in terms of its
% autocorrelation function.  Thus, 

% M4(i,j,k) = E { x(n) x(n+i) x(n+j) x(n+k) } 
%           = M2(i) * M2(j-k) + M2(j) * M2(k-i) + M2(k) * M2(i-j) 

% If x(n) is non-Gaussian - and most real world signals are non-Gaussian - 
% then, it is not completely characterized by its autocorrelation function.
% The higher-order moments of the process carry information which is not
% contained in the autocorrelation function. 

% Some examples of real-world signals that are not Gaussian are: 
%   speech  - Laplacian p.d.f.'s
%   radar   - Rayleigh distributions
%   seismic - generalized Gaussian distributions 
%   sonar   - mixture distributions 
%   bio-medicine, astronomy, image processing, oceanography, plasma physics ..

% Hit any key to continue
pause
clc
% So why should I use cumulants instead of moments ?  
% Several reasons: 
% 1. The cumulants of a Gaussian process are identically
%    zero for orders greater than two, e.g., 
% 
%        C3(i,j)   = 0 ,    C4(i,j,k) = 0,     etc 
% 
% 2. If z(n) = x(n) + y(n), where x(n) is statistically independent of x(n),
%    then, the cumulant of z(n) is the sum of the cumulants of x(n) and y(n)
% e.g., 
%     C2z(i)       = C2x(i)     + C2y(i) 
%
%     C3z(i,j)     = C3x(i,j)   + C3y(i,j) 
%
%     C4z(i,j,k)   = C4x(i,j,k) + C4y(i,j,k) 
%
% where C2z denotes the second-order cumulant of process z(n), etc.
% Note that the above property holds for second- and third-order moments
% of zero-mean processes, but not for higher-order moments. 

% Hit any key to continue
pause
clc
% 3. If z(n) = s(n) + g(n), with  `signal' s(n) statistically independent of 
%    `noise' g(n), and if g(n) is Gaussian, then, 
% 
%     C2z(i)     = C2s(i) + C2g(i)
%
%     C3z(i,j)   = C3s(i,j) 
%
%     C4z(i,j,k) = C4s(i,j,k) 
%                                .... etc. 
% Note that the effects of the Gaussian noise have been suppressed in 
% the higher-order cumulant domains! 
%
% Also C3z(i,j) = C3s(i,j) holds if the noise, g(n), is 
%      symmetrically distributed (e.g., Gaussian, Laplace, Uniform, etc.)
%
% What about the signal information? 

% Hit any key to continue
pause
clc

% 4. Suppose x(n) is an i.i.d. process.
%    Then, its cumulants are non-zero only at the origin, i.e., 
%
%    C2(n)   == 0, for n not 0
%    C3(m,n) == 0,  for (m,n) not equal to (0,0), etc. 
%
%    This property is useful in analyzing linear processes. 
%
%    The moment functions do not enjoy this nice property. 

% Hit any key to continue
pause
clc

% 5. Suppose x(n) is a linear non-Gaussian process, i.e, 
%
%         x(n) = sum {over i} h(i) u(n-i) 
%
% where u(n) is an i.i.d. non-Gaussian process, then, 
% the cumulants of x(n) can be expressed in terms of impulse response h(n),
% 
%    C2x(i)      =  C2u(0)     sum {over n} h(n) h(n+i) 
%    C3x(i,j)    =  C3u(0,0)   sum {over n} h(n) h(n+i) h(n+j)
%    C4x(i,j,k)  =  C4u(0,0,0) sum {over n} h(n) h(n+i) h(n+j) h(n+k)
%
% Again, this simple property does not hold for moments. 
% If C3u(0,0) is not zero, then under mild conditions, C3x(i,j) contains 
% sufficient information to recover both the magnitude and phase of H(f), 
% the Fourier transform of h(n).
% Recall that if we use C2x(i), we can only recover the magnitude; to recover 
% the phase, we must make an ASSUMPTION, such as minimum-phase or zero-phase. 
% Using C3(i,j) or C4(i,j,k), we no longer need to make the minimum-phase
% assumption.  The problem of estimating h(n), given only x(n),
% is called `blind deconvolution', and arises in a wide variety of
% situations, for example, in seismic signal processing, in 
% communications signal processing (channel effects),  
% image restoration (blurring), and so on. 

% Hit any key to continue
pause
clc

% 6. Non-linear processes: 
% Recall that from the Wold decomposition, any finite energy, wide-sense
% stationary, strictly random process can be modeled as a linear process 
% with respect to its second-order statistics: i.e., given x(n), 
% there exists a finite energy filter, h(n), and an uncorrelated (white)
% process, u(n), with finite variance, such that 
%
%    x(n) = sum {over i} h(i) u(n-i)
%
% In other words, second-order statistics do not yield any information
% about non-linear mechanisms that may have generated x(n).

% Hit any key to continue
pause
clc

% 6. Non-linear processes:  (continued)
% Non-linear phenomena are encountered frequently in real world situations.
% For example,  in fluid mechanics, plasma physics, ocean wave couplings,
%   EEG signals, otoacoustic emissions,  rotating machinery, etc. 
% The general theory of Volterra systems describes these phenomena.
% Let us concentrate on zero-memory power law mechanisms, such as 
% the effects of a squaring operation or a cubing operation. 
%
%   y(n) = x(n) ^ 2        z(n) = x(n) ^ 3 
%
% If x(n) consists of two harmonics at frequencies f1 and f2, then y(n) 
% contains harmonics at f1, f2, f1+f2 and f1-f2.   The distortion 
% harmonics at f1+f2, and f1-f2 are said to be frequency (phase) coupled
% to the original harmonics at f1 and f2.   This phenomenon is also called
% Quadratic Phase Coupling (QPC), and is encountered in plasma physics,
% oceanography, EEG signals, otoacoustic emissions, signals generated
% by rotating machinery, etc.  
%
% The presence of quadratic or cubic phase-coupling cannot be detected 
% using the autocorrelation.  But, they can be detected and quantified
% using the third- and fourth-order cumulants. 

% Hit any key to continue
pause
clc

% 6. Non-linear processes:  (continued)
% 
% Volterra series or representations are useful to represent the input
% output mapping in a non-linear system.
% The simplest non-linear system in this representation is the so-called
% second-order Volterra model,  described by 
% 
%  y(n) =   sum (over k)   h(k) x(n-k) 
%         + sum (over k,l) q(k,l) x(n-k)x(n-l)
% 
% The quadratic kernel is usually assumed to be symmetric.
% Note that a least-squares formulation to estimate the h(k)'s, and 
% the q(k,l)'s will involve second-, third-, and fourth-order moments.
% The QPC problem is a special case where the linear part, h(k), is 0,
% the quadratic part, q(k,l), is diagonal, and x(n) is a sum of harmonics. 

% Hit any key to continue
pause
clc

% 7.  Multi-channel signals: spatially correlated noise and cross-cumulants 

% So far we have considered scalar signals.  In various real-word scenarios,
% such as in sonar or seismic, the received signal is recorded at an
% array of sensors. The signal is corrupted by noise which may be 
% both spatially as well as temporally correlated.  
% Consider the Time Delay Estimation (TDE) problem: 
% here, we want to find the difference in the arrival time of a signal
% at two different sensors (this is the basic problem in the 
% Direction of Arrival (DOA) problem, as well, except, that multiple
% sensors are used). 
% If the signals are corrupted by spatially correlated noise, second-order 
% statistics are inadequate, since we cannot discriminate the delay due to
% the signal form that due to the noise.  
 
% Hit any key to continue
pause
clc 

% 7.  Multi-channel signals: spatially correlated noise and cross-cumulants 
%     ... continued

% Let      y(n,i) = x( n+d(i) ) + g(n,i),  i = 1,2,3 .... 
%
% where x(n) is the basic signal, d(i) is the delay with respect to some
% arbitrary reference point, and g(n,i) is spatially correlated Gaussian noise
% at time n and sensor i. 
% Assume without loss of generality, that the signals and noise are zero-mean,
% and that the signal is independent of the noise.   Then, 
%
%  E { y(n,i) y(n+m,j) y(n+p,k) } = C3x( m+d(j)-d(i) , p+d(k)-d(i) ) 
% 
% cum{ y(n,i), y(n+m,j), y(n+p,k), y(n+q,l) } 
%      = C4x( m+d(j)-d(i), p+d(k)-d(i), q+d(l)-d(i) )
%      
% cum{ y(n,i), y(n+m,i), y(n+p,i), y(n+q,i) }  = C4x( m, p, q )
% Note that the spatially correlated Gaussian noise has been suppressed. 

% The relative time-delays,  d(j)-d(i) , can be recovered from the
% auto- and cross-cumulants. 

% Hit any key to continue
pause
clc


%                 The bispectrum, and the bicoherence 
% 
% The BISPECTRUM is defined as the 2-D Fourier Transform of the third-order
% cumulant, 
%           B(f1,f2) = sum {over m,n} C3(m,n) exp( -j 2 pi (m f1 + n f2) )
%
% Note that the bispectrum is a function of two frequency variables. 
% Similarly, the trispectrum is defined as the 3-D Fourier transform of the
% fourth-order cumulant, and is a function of three frequency variables. 
% 
% If x(n) is a linear process, i.e.,    x(n) = sum {over i} h(i) u(n-i)
% where u(i) is a white process, then
%    Sx(f)      = C2u(0)    H(f)  conj ( H(f) ) 
%    Bx(f1,f2)  = C3u(0,0)  H(f1) H(f2)  conj ( H(f1+f2) ) 
% where Sx(f) denotes the power spectrum of process x(n). 
%
% The bicoherence is defined by 
%  Bic (f1,f2) =  B(f1,f2) / sqrt{ S(f1) S(f2) S(f1+f2) } 
% 
% Note that the absolute value of the bicoherence is a constant for 
% the linear process x(n). 

% Hit any key to continue
pause
clc

%            The cross bispectrum and the cross-bicoherence 
% 
% The cross-cumulant of three zero-mean random processes, x(n), y(n) and z(n)
% is defined by
%                 Cxyz(m,n) = E { x(i) y(i+m) z(i+n) } 
%
% The cross-bispectrum is the 2-D Fourier Transform
% 
%    Bxyz(f1,f2) = sum {over m,n} Cxyz(m,n) exp( -j 2 pi (m f1 + n f2) )
%
% The cross-bicoherence is defined by 
%  BICxyz (f1,f2) =  Bxyz(f1,f2) / sqrt{ Sx(f1) Sy(f2) Sz(f1+f2) } 
% 
% The cross-cumulant of four random processes, w(n), x(n), y(n), and z(n)
% is defined by
% 
%  Cwxyz(k,m,n) = cum(w(i), x(i+k), y(i+m), z(i+n)) 
%               = E{ w(i)x(i+k)y(i+m)z(i+n) } - E{w(i)x(i+k)} E{y(i+m)z(i+n)}
%                -E{w(i)y(i+m)} E{x(i+k)z(i+n)} -E{w(i)z(i+n)} E{y(i+m)x(i+k)}
%
% The fourth order cross-spectrum and tricoherence are similarly defined.
% 
% Hit any key to return to the main menu
pause
echo off
clc

