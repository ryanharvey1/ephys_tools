%DEMVOLT HOSA Toolbox Demo:  Identification of Second-Order Volterra Systems

echo off
% demos   nltick, nlpow  (nlgen)

% A. Swami, Nov 15, 1994
% Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.8 $

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

%       Identification of Second-Order Volterra Systems
%
% The HOSA Toolbox offers two routines for estimation of the parameters of a
% second-order Volterra system, given both the inputs and outputs:
%
% NLTICK - uses an algorithm due to Tick  (assumes Gaussian inputs)
% NLPOW  - uses an algorithm due to Powers (arbitrary inputs)
%
% The output of a second-order Volterra system is given by,
%  y(n) =   sum (over k)   h(k) x(n-k)
%         + sum (over k,l) q(k,l) x(n-k)x(n-l)
%
% Routine  NLGEN can be used to compute the output.
%
% If the input process, $x(n)$, is Gaussian, then the linear part { h(k) } can
% be determined from the cross-spectrum between x(n) and y(n).
% The quadratic part, { q(k,l) }, which is assumed to be symmetric can be
% estimated from the  cross bispectrum between processes x(n) and y(n),
% using HOSA routine BISPECDX.

%Hit any key to continue
pause

load nl1

% A white Gaussian process, x, was passed through a second-order Volterra
% system, to obtain the process y.   The columns of the 64 by 64 matrices
% x and y correspond to different realizations.
% The linear part { h(k) } is the IR of a low-pass MA(11) filter, with a
% nominal cutoff of 0.2 Hz:   h = fir1 (11, 0.4 )
% The quadratic part { q(k,l) } was chosen as q(k,l) = g(k)g(l), where
% g = fir1(11, 0.2).
% We used NLGEN to generate output y, from input x, and filters h and q.


echo off
h = fir1(11, 0.4);
g = fir1(11, 0.2);  q = g' * g;

clf,
nfft=128; w = [0:nfft/2]/nfft;
hf = fft(h,nfft);
subplot(221)
semilogy (w, abs(hf(1:nfft/2+1)) ),  title('linear TF'), grid on
q1 = zeros(32,32);
n = length(g);
q1(1:n,1:n) = q;
qf = fftshift( fft2(q1,nfft,nfft) );
subplot(222)
w2 = [-nfft/2:nfft/2-1]/nfft;
%contour(abs(qf), 6, w2, w2); title('quadratic TF')
contour(w2,w2,abs(qf), 6);  title('quadratic TF'), grid on
subplot(223), plot(h), title('linear part: IR'), grid on
subplot(224), contour(q1),  title('quadratic part: IR'), grid on
drawnow
set (gcf, 'Name', 'HOSA - True Volterra model')

echo on
%
% Hit any key to continue
pause

echo off

l_volt = str2mat( ...
         'NLTICK - Gaussian inputs ', ...
         'NLPOW  - Arbitrary inputs ');
c_volt = str2mat('d_nltick','d_nlpow');

choices ('HosatVoltDemo',' HOSA - Volterra Systems',l_volt, c_volt, 1);
echo off
return