%D_NLPOW  HOSA Demo: Identification of Second-Order Volterra Systems
%                      Powers' algorithm, for arbitrary inputs 
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

%  	NLPOW:     Powers' algorithm, for arbtirary inputs 
%
%   Routine NLPOW can  be used to estimate the parameters of a second-order
%   Volterra system; it does not make any assumptions about process x(n). 
%   We will use NLPOW on the data we just used.
%

%Hit any key to continue
pause 

    load nl1 
    clf
    [hp,qp] = nlpow(x,y,128);
    set (gcf, 'Name', 'HOSA NLPOW')
    drawnow

%   Hit any key to return to the main menu
   pause 
   echo off
clc 

