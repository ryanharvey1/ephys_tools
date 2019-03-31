%D_ORDER  HOSA Demo:  Linear Processes - ARMA model order determination.
echo off
% Demo of  arorder, maorder


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

%              Model order determination
%
% HOSA offers the following routines related to order determination
%
% ARORDER - estimates AR order based on the singular values of a cumulant or
%	    correlation matrix.
% MAORDER - estimates the MA order based on a hypothesis testing procedure
%	    applied to third-order cumulants.
%
% Hit any key to continue
pause
echo off
l_order = str2mat('AR order estimation ', ...
	          'MA order estimation ');
c_order = str2mat('d_arord', 'd_maord');
choices('HosatOrderDemo','HOSA - ARMA order determination ', ...
	  l_order, c_order, 1);

echo off
return