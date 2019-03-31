%EDA
% data in sp:   sp_txt: data description (string)
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
% 

if (exist('rlags') ~= 1) rlags = 30; end     % correlation lags:harmest
%       $Revision: 1.6 $
rlags = min(length(sp)-1, max(2,rlags) );  
if (exist('clags') ~= 1) clags = 25; end     % cumulant lags: bispeci
clags = min(length(sp)-1, max(2,clags) );  
if (exist('rfft')  ~= 1)  rfft =256; end     % fft length
if (exist('do_pcolor') ~= 1) do_pcolor = 1;  end 

 c1 = mean(sp);
 c2 = cumest(sp,2);
 c3 = cumest(sp,3) / c2^(3/2);
 c4 = cumest(sp,4) / c2^2;


% -------- Plot data and histogram 
 figure
 figno(length(figno)+1) = gcf; 
 subplot(211),  plot( 1:length(sp),sp), grid on
 subplot(212), hist(sp), 
 set(gcf,'Name',sp_txt)

 fprintf('Data and histogram plotted in figure window\n');
 fprintf(' ------  Summary stats \n'); 
 fprintf(' Mean                    %g\n',c1);
 fprintf(' Variance                %g\n',c2);
 fprintf(' Skewness (normalized)   %g\n', c3);
 fprintf(' Kurtosis (normalized)   %g\n\n', c4); 

 disp('Hit any key to continue')
 pause
 clc  
 echo on 

% -------- power spectra and harmonic models:   
 figure
 figno(length(figno)+1) = gcf; 
 [p2,a2,b2] = harmest(sp, rlags, 0, 'u', rfft, 2);    
 echo off 

 set(gcf,'Name',[sp_txt '- power spectra'])
 r = cplxpair(roots(a2));
 a = angle(r);
 yest = a*0;
 ind = find(a > 0) ; 
 yest(ind) = (2*pi) * ones(size(ind)) ./ angle(r(ind))  ;
 disp('Estimated cycles')
 disp(yest(ind))

 echo on 
% Hit any key to continue 
 pause 

% -------- gaussianity and linearity tests 
 glstat(sp, .51, length(sp) ) ;

% If the PFA is very small, the data are consistent with the hypothesis of
% non-zero bispectrum.  If, in addition, the estimated and theoretical values
% of 'R' are very different from one  another, then, the linearity hypothesis
% must be rejected. 

  
% -------- the bispectrum 
 figure
 figno(length(figno)+1) = gcf; 
 [Bspec, w] = bispeci(sp, clags);     




echo off 

 if (do_pcolor) 
    pcolor(w,w,abs(Bspec)), shading('interp')  % nice but slow 
 end
 set(gcf,'Name',[sp_txt '- bispectrum'])

disp('Hit any key to continue')
pause 
clc
echo on
