%EX_CASE: HOSA Toolbox Demo: Case studies
%
echo off

% A. Swami April 15, 1995
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

clear, 
echo on

% We will use some of the HOSA routines to analyze some real data:
%    sunspot data
%    canadian lynx data
%    speech data - a segment of ``laughter''

%Hit any key to continue
pause
echo off

l_case = str2mat('Sunspot data ','Canadian lynx data','Speech data');
c_case = str2mat('ex_suns', 'ex_lynx', 'ex_laff' );

choices('HosatCaseDemo','HOSA - Case studies', l_case, c_case,1);
echo off
clc
return