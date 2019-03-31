%D_RIVTR HOSA Demo: Adaptive AR parameter estimation - transversal form
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

%        Adaptive AR parameter estimation - transversal form
%
%    RIVTR obtains adaptive estimates of the AR parameters, using the 
%    transversal form of the recursive instrumental variable algorithm. 
%    The algorithm generates time-varying AR parameters and the 
%    final prediction error. 

% Hit any key to continue
pause 

% Now, we will use RIVTR to estimate the AR parameters. 
% 
load riv
[ar(:,1),fpe,wt] = rivtr(y,2);      % second-order cumulants
ar(:,2) = rivtr(y,3);               % third-order cumulants
ar(:,3) = rivtr(y,4);               % fourth-order cumulants
ar(:,4) = rivtr(zw,4);		    % fourth-order cumulants
ar(:,5) = rivtr(zc,4);              % fourth-order cumulants
ar(:,6) = rivtr(zc,2);              % second-order cumulants


% The true AR parameters were [1, -1.5, 0.8]. 
% The Steady-state AR estimates obtained using RIVTR are: 
disp(ar) 
% The first three columns correspond to estimates obtained from the noise-free
%     signal, and are based on the second, third, and fourth-order cumulants. 
% The fourth column corresponds to the signal corrupted by white Gaussian 
%     noise, and is based on the fourth-order cumulant. 
% The fifth and sixth columns correspond to the signal corrupted by colored 
%  Gaussian noise, and are based on the fourth-order cumulant, and the 
%  correlation, respectively.   In the last case, note that the results
%  are badly biased: a second-order model is inadequate. 

% Let us look at the temporal evolution of the weight vectors 
%  corresponding to the first test case, in order to see how fast the
%  estimates converge

% Hit any key to continue
pause 

clf
subplot(211)
plot(wt),title('weight vectors'), grid on
set(gcf,'Name','HOSA RIVTR')

% Hit any key to return to the previous menu 
pause 
echo off
clc
