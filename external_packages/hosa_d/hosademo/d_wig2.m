%D_WIG2  HOSA Demo of Wigner Spectrum  (wig2 and wig2c) 
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

%                      Wigner Spectrum 
%

% Now let us compute the conventional Wigner spectra:
% Hit any key to continue
pause

clf 
load wigdat 
subplot(221), wig2(s2,[ ],0);  title('WS - a')
subplot(222), wig2(s2);        title('WS - b')
subplot(223), wig2(s3);        title('WS - c')
subplot(224), wig2(s2+s3);     title('WS - d')
set(gcf, 'Name', 'HOSA WIG2')

% Note that the Wigner Spectrum of signals s2 and s3 are concentrated around 
% their nominal center frequencies. 
% 
% In panel (a) note the interaction between the energies at the positive
% and negative frequencies giving rise to the interference terms around D.C. 
% In panel (b) we use the analytic version of the signal so that the 
% negative frequency terms are suppressed;  this also suppresses the 
% distortion terms around D.C. 

% In panels (b) and (c) note that the Wigner spectrum is concentrated
% around the nominal center frequencies.

% In panel (d), we show the Wigner spectrum of the sum signal. 
% Note the cross-terms: some of these can be eliminated by using routine WIGC

% Now we will attempt to suppress the cross-terms by using the Choi-Williams
% filtering in routine WIG2C. 
% The amount of cross-term suppression is dictated by the filter parameter
% sigma: large values are equivalent to no filtering;  very small values will
% result in signal distortion. 

% Hit any key to continue
pause

subplot(221)
wig2c(s4); title('wig2c - a')

% The WIG2C output is plotted in panel (a).
% Hit any key to return to previous menu
pause
echo off
clc 
return

