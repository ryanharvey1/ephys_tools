%DEMGLTST HOSA Toolbox Demo:  Tests for Gaussianity and Linearity 

echo off 

% A. Swami April 15, 1993
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

%             TESTING for GAUSSIANITY and LINEARITY 
%
% How do we know whether `real' data are non-Gaussian or non-linear? 
% 
% The basic idea is that if the signal is Gaussian, its third (fourth ...)
% order cumulants must be identically zero.  In practice, sample estimates
% of cumulants will not be exactly zero: so, we need a test
% to determine whether or not estimated quantities are significantly 
% different from zero in a statistical sense. 
% For a linear non-Gaussian process, we saw that the absolute value
% of the bicoherence is a constant.   Again, sample estimates of the 
% bicoherence will not be constant, and we need a test to
% determine whether the non-constancy is statistically significant. 

% The HOSA Toolbox offers the routine GLSTAT to test whether a given signal
% is non-Gaussian, and if so, whether it is also linear. 

% Hit any key to continue
pause
clc 

% Decision Statistics for Linearity and Gaussianity Tests 
%
load gldat 
% In this routine, the bispectrum of the process is estimated and 
% smoothed;   tests are then conducted to see whether the bispectral
% values are significantly different from zero.   The basic idea is that
% estimates of the bispectrum are asymptotically complex normal;  hence, the
% energy in the bispectrum is chi-squared distributed;  the number of 
% degrees of freedom depend upon the FFT length and the smoothing window. 
%
% In the Gaussianity test, 
% The null hypothesis is that the data have zero bispectrum ("Gaussian") 
% The computed probability of false alarm (PFA) value is the probability that 
% the value of the chi-squared r.v. with the indicated degrees of freedom will
% exceed the computed test statistic. 
% The PFA value indicates the false alarm probability in accepting the 
% alternate hypothesis, that the data have non-zero bispectrum. 
% Usually, the null hypothesis is accepted if PFA is greater than 0.05 
%    (i.e., it is risky to accept the alternate hypothesis). 

% Hit any key to continue
pause

% In the Linearity test, the inter-quartile range of the estimated
% bicoherence is computed;  a quantity, 'lambda' proportional to the
% mean value of the bicoherence is also computed;  the theoretical
% inter-quartile range of a chi-squared r.v. with two degrees of freedom
% and non-centrality parameter 'lambda' is then computed. 
% The linearity hypothesis should be rejected if the estimated and
% theoretical inter-quartile ranges are very different from one another.

% Hit any key to continue
pause 

% We will use a smoothing parameter (cparm) value of 0.51 and an FFT length
% of 256 in the following examples.  Each of the sequences to be tested 
% has 512 samples. 

      cparm = 0.51;  nfft = 256; 

% hit any key to continue 
pause 
clc

% We will apply the test to an i.i.d. Gaussian sequence, g. 

    [sg,sl] = glstat(g, cparm, nfft); 

% Since the PFA is high, we accept the null (Gaussian) hypothesis; 
% the linearity test is also based on the bispectrum; if the bispectrum is
% zero, the bicoherence will be a constant, equal to zero, and we cannot
% conclude, based on the bispectrum, whether or not the data are linear; 
% hence, the linearity test is meaningless in this case. 

% hit any key to continue  ............... 
pause 

% We will apply the test to an i.i.d. sequence, u, with uniform p.d.f.

    [sg,sl] = glstat(u, cparm, nfft); 

% Since the PFA is high, we accept the null (zero bispectrum) hypothesis;
% the linearity test is also based on the bispectrum; if the bispectrum is
% zero, the bicoherence will be a constant, equal to zero, and we cannot
% conclude, based on the bispectrum, whether or not the data are linear; 
% hence, the linearity test is meaningless in this case. 

% hit any key to continue  ............... 
pause 

% We will apply the test to an i.i.d. exponential sequence, e 

    [sg,sl] = glstat(e, cparm, nfft); 

% Since the PFA is very small, we accept the alternate hypothesis, 
%     i.e., the data are accepted as being non-Gaussian. 
% The linearity test is meaningful in this case. 
% The estimated and theoretical inter-quartile ranges are close to each other.
% Hence, we accept the linearity test as well. 

% hit any key to continue  ............... 
pause 

%  Sequence x was obtained by passing e through a linear filter. 
%  Since sequence  e  was accepted as non-Gaussian and linear, 
%     we expect x to be accepted as non-Gaussian and linear as well. 
%  Let us apply the tests to the sequence e. 

    [sg,sl] = glstat(x, cparm, nfft); 

% Since the PFA is very small, we accept the alternate hypothesis
%     i.e., the data are accepted as being non-Gaussian. 
% The linearity test is meaningful in this case. 
% The estimated and theoretical inter-quartile ranges are close to each other.
% Hence, we accept the linearity test as well. 

% hit any key to continue  ............... 
pause 

% Sequence z was obtained by passing the sequence x through a non-linearity, 
%   z = x.^3  
% Let us apply the tests to z 

     [sg,sl] = glstat(z, cparm, nfft); 

% Since the PFA is small, we accept the non-Gaussian hypothesis 
% Since the estimated and theoretical inter-quartile ranges are very 
% different, we cannot accept the linearity hypothesis. 

% hit any key to continue  ............... 
pause 

% We will apply the test to an i.i.d. Laplacian sequence, l  

    [sg,sl] = glstat(l, cparm, nfft); 

% Since the PFA is high, we accept the null (zero bispectrum) hypothesis;
% the linearity test is also based on the bispectrum; if the bispectrum is
% zero, the bicoherence will be a constant, equal to zero, and we cannot
% conclude, based on the bispectrum, whether or not the data are linear; 
% hence, the linearity test is meaningless in this case. 

% Hit any key to return to the main menu .... 
pause 
echo off
clc
