%D_RIVDL HOSA Demo: Adaptive AR parameter estimation - Lattice form
%       
echo off
% A. Swami Oct 18, 1997.
% Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.7 $

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

% 	   Adaptive AR parameter estimation - Lattice form
%
%    RIVDL obtains adaptive estimates of the AR parameters using the
%    double lattice form of the recursive instrumental variable algorithm. 
%    The estimates may be based on second, third, or fourth-order cumulants. 
%    The algorithm generates time-varying AR parameters, as well as the 
%    forward and backward reflection coefficients (from each stage, as 
%    functions of time) from the upper lattice. 
%

load riv

% We will use RIVDL to estimate the AR parameters. 
%    this routine tends to be a bit slow  ..... 

% Hit any key to continue 
pause 

[ar(:,1),fref,bref] = rivdl(y,2);     % second-order cumulants
ar(:,2) = rivdl(y,3);                 % third-order  cumulants 
ar(:,3) = rivdl(y,4);                 % fourth-order cumulants
ar(:,4) = rivdl(zw,4);                % fourth-order cumulants, white noise
ar(:,5) = rivdl(zc,4);                % fourth-order cumulants, colored noise
ar(:,6) = rivdl(zc,2);                % second-order cumulants, colored noise

% The true AR parameters were [1, -1.5, 0.8] 
% The steady-state AR parameters estimated by RIVDL are: 
disp(ar) 
% The first three columns correspond to estimates obtained from the noise-free
%   signal, and are based on the second, third, and fourth-order cumulants. 
% The fourth column corresponds to the signal corrupted by white Gaussian 
%   noise, and is based on the fourth-order cumulant. 
% The fifth and sixth columns correspond to the signal corrupted by colored 
%  Gaussian noise, and are based on the fourth-order cumulant, and the 
%  correlation, respectively.   In the last case, note that the results
%  are badly biased: a second-order model is inadequate. 

%
% Let us look at the temporal evolution of the reflection coefficients, 
%  corresponding to the first test case, in order to see how fast the
%  estimates converge

% Hit any key to continue
pause 


clf, 
subplot(211)
plot(fref),title('Forward  reflection coefficients'), grid on
subplot(212)
plot(bref),title('Backward reflection coefficients'), grid on 
set(gcf,'Name','HOSA RIVDL')

% Hit any key to return to the previous menu 
pause 
echo off
clc 
