%D_WIG3  HOSA Demo of Wigner Bispectrum  (wig3 and wig3c) 
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

%                    Wigner Bispectrum
%
% Now let us look at the Wigner bispectrum: 
% Routine WIG3 computes the f1=f2 slice of the Wigner bispectrum W(t,f1,f2) 

% Hit any key to continue
pause

load wigdat 
clf
subplot(221), wig3(s2,[ ],0);    title('WB - a')
subplot(222), wig3(s2);          title('WB - b')
subplot(223), wig3(s3);          title('WB - c')
subplot(224), wig3(s2+s3);       title('WB - d')
set(gcf, 'Name', 'HOSA WIG3')


% In panel (a) we did not use the analytic form of the signal; 
% hence, interaction between the energies at the positive and negative
% frequencies gives rise to the interference terms around D.C. 
% These are suppressed by using the analytic form in panel (1,2). 

% In panels (b) and (c) note that the Wigner bispectra are 
% concentrated around their nominal center frequencies. 

% In panel (d) we show the Wigner bispectrum of the signal s2+s3; note the 
% cross-terms; some of these can be eliminated by using routine WIG3C. 

% Now we will attempt to suppress the cross-terms by using the Choi-Williams
% filtering in routine WIG3C. 
% The amount of cross-term suppression is dictated by the filter parameter
% sigma: large values are equivalent to no filtering;  very small values will
% result in signal distortion. 

% Hit any key to continue
pause

subplot(221) 
wig3c(s4,256,0.2,1);  title('wig3c - a')

% The WIG3C output is plotted in panel (a)
% It should be noted that unlike WIG2C and WIG4C, WIG3C may distort the
% signal badly.

% Hit any key to continue
pause
echo off
clc 
return


