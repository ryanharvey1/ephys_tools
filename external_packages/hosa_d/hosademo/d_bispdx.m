%D_BISPDX  HOSA Demo: Estimating the Cross-Bispectrum (bispecdx)
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
%          BISPECDX - cross-bispectrum estimation  (direct method) 
%
% Here, the given time-series, x(n), y(n), z(n), are segmented into several, 
% possibly overlapping records,  xi(n), yi(n), zi(n), i=1,...,K, n=1,...,M. 
% The Fourier transforms of the records are computed, Xi(f), Yi(f), Zi(f).
% The sample triple product is computed, Xi(f1) Yi(f2) conj(Zi(f1+f2)), 
% and averaged across the set of records,
%   1/MK sum {over i} Xi(f1) Yi(f2) conj( Zi(f1+f2) ) 
% A frequency-domain smoothing window is then applied to obtain the
% final estimate. 
% 
%
% Let x(n) be a white Gaussian process, and let y(n) = x(n) .^2
% It is easy to show that the cross-bispectrum of (x,x,y)
% should be a constant, equal to 2 ( std(x) ) ^4 . 
%
% Hit any key to continue 
pause

clf
x = randn(64,64);
y = x.^2;
dbic = bispecdx(x,x,y,128,5); 
set (gcf, 'Name','HOSA - BISPECDX')

% The apparent structure along the axes is an artifact due to the 
% removal of the mean:  consistent estimates along the axes and the 
% anti-diagonal can be obtained only by sufficient smoothing; a smoothing
% window of size 5 is inadequate.  
% (The lines f1=0, f2=0, and f1+f2=0 are called the principal sub-manifolds). 
%
% Hit any key to return to main menu
pause 
echo off
clc

