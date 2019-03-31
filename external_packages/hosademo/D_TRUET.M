%D_TRUET  HOSA Demo for computing true trispectrum (trispect)
echo off

% A. Swami Jan 20, 1995 
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

% TRISPECT:  computes 1-D slices of the theoretical trispectrum of an 
%		ARMA model
%
% We will compute several 1-D slices of the trispectrum of the ARMA model 
%         with AR = [1, -1.5, 0.8] and MA = [1,-2] 
% and then, we will display a `movie'. 

%Hit any key to continue
pause 

clf 
set(gcf,'Name','HOSA trispect') 
ma = [1 -2];  ar = [1 -1.5 0.8];  nfft = 64; 
     n=5; M = moviein(2*n+1); 
     for k=-n:n
         trispect(ma,ar,nfft,k/(2*n));
         M(:,k+n+1) = getframe;
     end
     movie (M)
     clear M

% Hit any key to continue
pause 
echo off
clc
return

