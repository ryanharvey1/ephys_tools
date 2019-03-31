%EX_LAFF: Analysis of laughter data
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

 load laff.dat 		       % LAUGHTER data
 y = laff;   Fs = 8192;

echo off
%--- summary stats:

x=(y-mean(y))/std(y);
figno(1) = gcf; 
subplot(211), plot((1:1400)*1000/Fs, x), grid on 
xlabel('time in ms')
title('laughter data')
subplot(212), hist(x,20),

 fprintf('\n Data and histogram plotted in figure window\n');
 disp(' The marginal distribution does not appear to be symmetric')
 fprintf('\n ------  Summary stats \n');
 fprintf(' Mean                    %g\n',mean(y));
 fprintf(' Variance                %g\n',std(y).^2 );
 fprintf(' Skewness (normalized)   %g\n', mean(x.^3) );
 fprintf(' Kurtosis (normalized)   %g\n\n', mean(x.^4) - 3);

disp('Hit any key to continue')
pause

% The Spectrogram:
disp('We will now look at the spectrogram')
      figure
     figno(length(figno)+1) = gcf; 
  [spx, f,t] = specgram (x, 512, Fs, hamming(256), 240);
  contour(t*1000,f,abs(spx),8), grid on 
  xlabel('time in ms'), ylabel('frequency in Hz'),
  set(gcf,'Name','speech spectrogram')

echo on

% Note the dominant frequency tracks around 550, 1100 and 1600 Hz')
%     additional fragments are around 1800 Hz and 2100 Hz ')
% The data appear to be ``harmonic'' in nature.')

% Hit any key to continue
pause

% We can use parametric techniques to estimate the frequencies
% note that the segment is ``stationary''

% ``Power spectra'':
% The singular value plot indicate an order of p=12
     figure, [pxx,a2_1,a2_2] = harmest(x,25,0,'biased',512,2);
     figno(length(figno)+1) = gcf; 
     set(gcf,'Name','speech: power spectra')

% Hit any key to continue
     pause

% ``Fourth-order cumulant spectra''
% The singular value plot will indicate an order of p = 8
% whereas that based on the correlation indicated order p=12.
% A possible explanation is that some of the harmonics are well-modeled as
% narrow-band Gaussian, and some as narrow-band non-Gaussian.
%
    figure, [pxx,a4_1,a4_2] = harmest(x,25,0,'biased',512,4);
     figno(length(figno)+1) = gcf; 
     set(gcf,'Name','speech: cumulant spectra')

% Hit any key to continue
   pause

% We can test whether the bispectrum of the data are statistically
% non-zero by using GLSTAT

     glstat(x,.51,256);

% Since the Pfa is close to 0, we can be virtually certain that the
% data have non-zero bispectrum and hence are non-Gaussian.
% The  estimated and theoretical R values are not close to each other
% indicating that the data are not linear.

% Hit any key to continue
pause

% Both the spectra and cumulant spectra indicate the presence of harmonics
% In order to check for quadratic frequency coupling, we can use
% the spectrum, pick peaks, and check whether f1+f2=f3
% Or we can compute the non-parametric bispectrum, via BISPECD or BISPECI
% or we can use the parametric method in QPCTOR.

disp('[Bs,w]=bispecd( x, 256, 0, 100, 0); ')
echo off
load laff_b
figure,
load laff_b
contour(w,w,abs(Bs),4), grid on 
title('Bispectrum estimated via the direct (FFT) method')
xlabel('f1'), ylabel('f2')
set(gcf,'Name','speech: bispectrum')
figno(length(figno)+1) = gcf; 
echo on

% Note that only the first quadrant is displayed (sufficient)
% Since the data consists of several harmonics, approximately given by
% f_k = k f_0, the bispectrum consists of a mess of impulses.
% Since the underlying model is of a single fundamental, and its harmonics,
% it suffices to use the diagonal of the bispectrum:

% Hit any key to continue
pause

figure,
d=diag(Bs);
plot(w,abs(d)), grid on 
figno(length(figno)+1) = gcf; 
set(gcf,'Name','speech: diagonal slice of bispectrum')

% Hit any key to continue
pause

% We can locate the peaks by using PICKPEAK

[loc,val]=pickpeak(d,3,5);
disp( w(loc)'  )

% Hit any key to continue
pause

% We can use the parametric method QPCTOR to confirm the quadratic frequency
% coupling.
figure,
ar=qpctor(x,25,10,256,100,0,'biased');
set(gcf,'Name','speech: QPC detection')
figno(length(figno)+1) = gcf; 

% The dominant peak is at (0.0625, 0.0625) indicating the presence of
% a second harmonic at 0.1250 (approx 1024 Hz)
% Also, the contour plot shows peaks at (0.0625, k*0.0625),
% indicating `self-coupling' of higher-orders as well.

%  The diagonal slice of the k-th order moment spectrum is  defined by
%
%  S_k(f) := E (  [ X(f) ]^k [ X(kf) ]^*   )
%

% We can also check for self-coupling by looking at the  magnitude of the
% `frequency-scaled' product
%
%  S(k; f) =     X (f) *  X (kf)
%
% which will show an impulse at  f=f_0, if the data contains harmonics
%   with frequencies f_0 and kf_0

% Hit any key to continue
pause

figure 
echo off
kfft = 2048;
Xf = fft(x,kfft) ;
Xf = Xf(1:kfft/2);
if (exist('lf') ~= 1) lf = 5; end
for k=1:4
 L=fix(kfft/(2*k));
 S = Xf(1:L).^k  .* conj( Xf(1:k:k*L) ) ;
 S = abs ( filter(ones(lf,1), 1, S) );
 eval(['S' int2str(k) ' = S;']);
 eval(['subplot(22', int2str(k), ')'])
 w = [1:(L-1)] / kfft * Fs;
 semilogy(w, S(2:L) ),    %avoid DC
 a= axis;  axis([0, round(40/k)*100, a(3), a(4) ]) ;
 ylabel(int2str(k))
 txt1 = [ 'CX(f) .* X(' int2str(k) 'f)' ];
 eval(['title (txt1)'])
 grid on 
end
figno(length(figno)+1) = gcf; 
disp('These diagonal slices indicate self-coupling of various orders')

echo on
% 
%  End of laughter demo
%  Hit any key to clear all plots, and return to previous menu
%
pause

echo off 
for k=1:length(figno)
  close(figno(k))
end
clear figno
