%D_BIC  HOSA Demo: Estimating the Bicoherence (bicoher)
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
%        BICOHER - Auto-bicoherence estimation  (direct method) 
%
% The bicoherence is defined by 
%
%  Bic (f1,f2) =  B(f1,f2) / sqrt{ S(f1) S(f2) S(f1+f2) } 
%
% where B(f1,f2) is the bispectrum and S(f) is the power spectrum. 
% For a symmetrically distributed process, the bicoherence is zero. 
% For a linear process, the bicoherence is constant.
% Non-constant values of the bicoherence  indicate non-linear phenomena.
% 
% Hit any key to continue
pause

% We will compute the bicoherence of the output of a second-order Volterra
% system, 
%   y(n) = sum_{k=0 to 11} h(k) x(n-k) 
%        + sum_{k=0 to 11} q(k,l) x(n-k)x(n-l)
% The filters h and q were generated via,
%   h = fir1(11,0.4); g = fir1(11,0.2);  q = g(:) * g(:)'
% A white Gaussian process, generated via RPIID, was chosen for x;
%   and the outputs, y, were computed via NLGEN.
% 
% We will use the default FFT length of 128,  and the default Hanning window;
% since we have multiple realizations, the segment size will be set to the row
% dimension of matrix zmat, and the overlap parameter will be set to 0.

% Hit any key to continue
pause

clf
load nl1
dbic = bicoher(y); 
set (gcf, 'Name','HOSA - BICOHER')


% The bicoherence is not constant, indicating non-linear phenomena. 
% The contour plot shows sharp peaks along the diagonal (f1,f2), which is
% indicative of process y(n) being generated via a squaring action;
% the non-sharpness of the peaks, as well as the presence of structure
% around the origin, indicates that the non-linear relationship is not
% purely of the form y = x.^2.   
% 

% Hit any key to return to main menu
pause 
echo off
clc

