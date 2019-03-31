function [yo, fo, to]=mtchg(varargin);

% Multitaper time-frequency coherences-gram
%
% This basically does the same thing as mtcsg, but scales
% the cross-spectra to become coherences

[x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers] = mtparam(varargin);
[y fo to] = mtcsg(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers);


nCh1 = size(y,3);
nCh2 = size(y,4);

yo = zeros(size(y));

% main loop
for Ch1 = 1:nCh1
	for Ch2 = 1:nCh2
		
		if (Ch1 == Ch2)
			% for diagonal elements (i.e. power spectra) leave unchanged
			yo(:,:,Ch1, Ch2) = y(:,:,Ch1, Ch2);
		else
			% for off-diagonal elements, scale
			yo(:,:,Ch1, Ch2) = abs(y(:,:,Ch1, Ch2).^2) ...
				./ (y(:,:,Ch1,Ch1) .* y(:,:,Ch2,Ch2));
		end
	end
end
			
% plot stuff if required

if (nargout<1)
	ImageMatrix(to,fo,yo);
end;

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu