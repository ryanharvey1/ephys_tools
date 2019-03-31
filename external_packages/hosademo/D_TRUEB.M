%D_TRUEB  HOSA Demo for computing true bispectrum (bispect)
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

% BISPECT:  computes the theoretical bispectrum of an ARMA model
%
% We will compute the bispectrum of the ARMA model 
%         with AR = [1, -1.5, 0.8] and MA = [1,-2] 

%Hit any key to continue
pause

clf 
bs = bispect ( [1 -2],  [1 -1.5 0.8], 256); 
title('True bispectrum: ARMA(2,1)')
set(gcf,'Name','HOSA BISPECT')

% Notice the  various symmetries in the bispectrum
% Note that this should be identical with the FT of the true cumulants
%    computed via CUMTRUE. 

%Hit any key to continue
pause
echo off
clc 

return
