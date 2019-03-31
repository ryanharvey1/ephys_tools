%D_WIG4  HOSA Demo of Wigner Trispectrum  (wig4 and wig4c) 
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

%		   Wigner Trispectrum
%
% Routine WIG4 computes the diagonal slice of the Wigner Trispectrum. 
% Let us look at some Wigner trispectra

% Hit any key to continue
pause

clf,  
load wigdat
subplot(221), wig4(s2,[ ],0);   title('WT - a')
subplot(222), wig4(s2); 	title('WT - b')
subplot(223), wig4(s3);		title('WT - c')
subplot(224), wig4(s2+s3); 	title('WT - d')
set(gcf, 'Name', 'HOSA WIG4')

% In panel (a), we used the real form of the signal; this leads to 
% interaction between the energies at the positive and negative frequencies,
% giving rise to the interference terms around D.C.   
%
% In panel (b) we use the complex form of the signal; this suppresses
% the negative frequencies and thus the distortion terms around D.C.
%
% In panels (b) and (c) note that the Wigner trispectrum is concentrated
% around the nominal center frequencies.
%
% Panel (d) shows the WT of the sum signal: note the cross-terms; some of
% these can be suppressed by routine WIG4C, as we shall see next. 

% Routine WIG4C attempts to suppress some of the cross-terms in the Wigner
% trispectrum, by applying a generalized form of the Choi-Williams filter.
% As in WIG4, the diagonal slice is computed.  This filtered diagonal slice
% has been called the ``Sliced Reduced-Interference Wigner Trispectrum''. 
% 

% Hit any key to continue
pause 

subplot(221)
wig4c(s4);  title('wig4c   - a')

% The WIG4C output is plotted in panel (a). 
% Note the suppression in cross-terms: 
% Hit any key to return to main menu 
pause 
echo off
clc
return

