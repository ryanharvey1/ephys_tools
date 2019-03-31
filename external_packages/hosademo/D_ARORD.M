%D_ARORD  HOSA Demo:  Linear Processes - AR model order determination.
echo off
% Demo of  arorder


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
%              ARORDER:  AR order determination 
%
% We will use ARORDER to determine the AR order of an ARMA(2,1) process
% which has been corrupted by additive white Gaussian noise, SNR = 20 dB. 
% We will use use both third and fourth-order cumulants: 
% ARORDER will choose the optimal order based on the difference in
% singular values; alternatively, we can inspect the singular value plot
% and choose it ourselves.  The plot of differences in singular values is
% expected to show a sharp peak at the true order. 

%Hit any key to continue
pause 

load arma1
clf
subplot(211), p3 = arorder(y,3);
subplot(212), p4 = arorder(y,4); 
set(gcf,'Name','HOSA ARORDER')

% Hit any key to continue
pause
echo off 
clc
