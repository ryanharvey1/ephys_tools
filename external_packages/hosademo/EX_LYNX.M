%EX_LYNX: Analysis of Canadian lynx data 
%
echo off

% A. Swami April 15, 1995 
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

clear, 
echo on


 load lynx.dat               % LYNX data 
 sp = lynx(:,2);

 figno(1) = gcf; 
 clc
 subplot(211),  plot( 1:length(sp),sp), grid on 
 subplot(212), hist(sp), 
 set(gcf,'Name','lynx data') 

%Hit any key to continue
pause  

% These are the  Canadian lynx data, counts for years 1821-1934
%   these data are all positive valued 
%   differencing the data is useful in such cases 
% The histogram (univariate) is almost symmetric, but the data have 
% non-zero bispectrum. 

 sp = diff(sp);
 sp_txt = 'lynx data (differenced)';
 ex_eda 

echo on
%  End of lynx-data demo
%  Hit any key to clear all plots, and return to previous menu
pause

echo off
for k=1:length(figno)
    close(figno(k))      % figno has figure numbers returned by ex_eda
end   
clear figno

return
