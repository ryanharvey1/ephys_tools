function [COH,freqVec, P]=signals2coh(S1,S2,Fs,nfft)
% function [COH,FREQVEC]=signals2coh(S1,S2,FS,NFFT)
% 
% Calculates the coherence between signals S1 and S1. A Hanning
% window of length nfft with overlap nfft/2 is applied. Signals are
% detrended prior to fft. 
%
% S1      : signal 1
% S2      : signal 2
% FS      : Sampling frequency
% NFFT    : number of fft points
% COH     : coherence
% FREQVEC : frequency vector
%
% See also: TRACES2PLF FIFF2COH
%
%------------------------------------------------------------------------
% Ole Jensen, Brain Resarch Unit, Low Temperature Laboratory,
% Helsinki University of Technology, 02015 HUT, Finland,
% Report bugs to ojensen@neuro.hut.fi
%------------------------------------------------------------------------

%    Copyright (C) 2000 by Ole Jensen 
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You can find a copy of the GNU General Public License
%    along with this package (4DToolbox); if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


w = hanning(nfft)';
ScaleFac = nfft/sum(w);

idxStart = 1;
idxEnd   = nfft;
k = 0;
while idxStart <= length(S1)
    k = k + 1;
    if idxEnd > length(S1)
        idxEnd = length(S1);
    end
    S1s = S1(idxStart:idxEnd);
    S2s = S2(idxStart:idxEnd);
    if length(S1s) < nfft
        S1s = [S1s zeros(1,nfft - length(S1s))];
        S2s = [S2s zeros(1,nfft - length(S2s))];

    end
    S1s = detrend(S1s);
    S2s = abs(detrend(S2s));
    PS1 = fft(w.*S1s,nfft);
    PS2 = fft(w.*S2s,nfft);
    PS1S2 = PS1.*conj(PS2);
    if exist('PS1sum')   
        PS1sum = PS1sum + abs(PS1).^2;
        PS2sum = PS2sum + abs(PS2).^2;
        PS1S2sum = PS1S2sum + PS1S2;        
    else
        PS1sum = abs(PS1).^2;
        PS2sum = abs(PS2).^2;
        PS1S2sum = PS1S2;        
    end
    idxStart = idxStart + floor(nfft/2);
    idxEnd   = idxEnd + floor(nfft/2);
end
P = (abs(PS1S2sum).^2)./(PS1sum.*PS2sum);
COH = P(1:floor(nfft/2));

freqVec = Fs*(0:floor(nfft/2)-1)/nfft;  

 



