%D_BICX  HOSA Demo: Estimating the Cross-Bicoherence (bicoherx)
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
%          BICOHERX - cross-bicoherence estimation  (direct method) 
%
% The cross-bicoherence between three processes, x,y and z, is defined by 
%
%  BICxyz (f1,f2) =  Bxyz(f1,f2) / sqrt{ Sxx(f1) Syy(f2) Szz(f1+f2) } 
%
% where Bxyz is the cross-bispectrum, and Sxx, Syy, Szz are the power spectra 
%
% The cross-bicoherence  will be zero if the three processes are jointly
%  symmetric distributed.
% A constant value of the cross-bicoherence is indicative of joint linear 
% relationships; 
% a non-constant value is indicative of non-linear relationships. 
% 
% Hit any key to continue
pause

%
% We will compute the cross-bicoherence between the inputs and outputs of
% a second-order Volterra system. 
% The data for the demo consist of the inputs, x, and outputs,y, of a 
% second-order Voltera system,
%   y(n) = sum_{k=0 to 11} h(k) x(n-k) 
%        + sum_{k=0 to 11} q(k,l) x(n-k)x(n-l)
% The filters h and q were generated via,
%   h = fir1(11,0.4); g = fir1(11,0.2);  q = g(:) * g(:)'
% A white Gaussian process, generated via RPIID, was chosen for x;
% and the outputs, y, were computed via NLGEN.

% Hit any key to continue
pause

clf
load nl1
dbic = bicoherx (x,x,y); 
set (gcf, 'Name','HOSA - BICOHERX')

% The cross-bicoherence is significantly non-zero, and non-constant, 
% indicating a non-linear relationship between  processes x and  y. 

% Hit any key to return to main menu
pause 
echo off
clc

