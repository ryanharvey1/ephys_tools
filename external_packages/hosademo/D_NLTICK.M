%D_NLTICK  HOSA Demo: Identification of Second-Order Volterra Systems
%                      Tick's algorithm, assuming Gaussian inputs
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

%      NLTICK:   Tick's algorithm, assuming Gaussian inputs
%
%   NLTICK assumes that the system input is Gaussian;  it computes the
%   cross-spectrum and cross-bispectrum of the system inputs and outputs to
%   estimate the linear and quadratic transfer functions. 
%   We will now use NLTICK to estimate the parameters, with an FFT length of
%   128; since we have multiple realizations, we will not use any windowing. 

%Hit any key to continue
pause 

    load nl1 
    clf
    [ht,qt] = nltick(x,y,128,1); 
    set (gcf, 'Name', 'HOSA NLTICK')
    drawnow

%   Hit any key to continue
    pause
echo off
clc 
return


