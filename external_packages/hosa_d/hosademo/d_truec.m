%D_TRUEC  HOSA Demo for computing true cumulants
%	demo of cumtrue

echo off

% A. Swami Jan 20, 1995
% Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.8 $

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
%	  Computing true (theoretical) cumulants
%
% CUMTRUE: computes the theoretical cumulants of ARMA models
%         second-, third- or fourth-order cumulants may be estimated
% We will compute the third-order cumulants of an ARMA(2,1) model
%         with AR = [1, -1.5, 0.8] and MA = [1,-2]
%

% Hit any key to continue
pause
cmat  = cumtrue([1,-2], [1,-1.5,0.8],3,25);

clf
subplot(121)
mesh(cmat),title('True cumulant: ARMA(2,1)')
subplot(122)
contour(-25:25,-25:25,cmat,8),title('True cumulant: ARMA(2,1)'), grid on
%contour(cmat,8,-25:25,-25:25),title('True cumulant: ARMA(2,1)'),
set(gcf,'Name','HOSA CUMTRUE')

% Hit any key to continue
pause
echo off
clc