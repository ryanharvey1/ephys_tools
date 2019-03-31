%D_MA  HOSA Demo:  Linear Processes - Parametric (MA) model estimation
%                        Blind deconvolution. 
echo off
% Demo of  maest 

% A. Swami Jan 15, 1995
% Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.5 $

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

% 
%           MA parameter estimation methods
%  
% MAEST estimates the MA parameters using an algorithm due to Giannakis-Mendel
% and a modification by Tugnait.    The algorithm uses both autocorrelation
% and cumulants, and assumes that the additive noise is white. 
%
% The test synthetic "x" is an MA(3) synthetic with true MA parameters, 
% [1, 0.9, 0.385, -0.771];  input was i.i.d. exponential; white Gaussian
% noise was added to obtain an SNR of 20 dB. 
%
%Hit any key to continue
pause 
% MA estimates based on third-order cumulants: 
load ma1

bvec = maest(y,3,3,128); 
disp(bvec')

% MA estimates based on fourth-order cumulants:  
bvec = maest(y,3,4,128); 
disp(bvec')

% Hit any key to continue
pause 
echo off 
clc

