function [yo, fo, to]=mtchgnorm(varargin);
% Multitaper time-frequency normalized coherences-gram
%
% This takes the output of mtcohere and applies transformations to
% try and make each variable normal.
%
% It does the following:
%
% For power spectra: it divides by the mean and multiplies by 2*nTapers.  
% You should then have a chi2 with 2*nTapers degrees of freedom
% which is converted to approximately normal by taking a square root
%
% For coherences it uses Fisher's z-transform

[x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers] = mtparam(varargin);
[y fo to] = mtchg(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers);

nFreqBins = size(y,1);
nTimeBins = size(y,2);
nCh1 = size(y,3);
nCh2 = size(y,4);

n = 2*nTapers; %degrees of freedom

yo = zeros(size(y));

% main loop
for Ch1 = 1:nCh1
	for Ch2 = 1:nCh2
		
		if (Ch1 == Ch2)
			% for diagonal elements (i.e. power spectra) normalize by 
			% mean power in each frequency band
			ynorm = y(:,:,Ch1, Ch2);
			ynorm = n*ynorm ./ repmat(mean(ynorm,2), 1, nTimeBins);
			% and then do chi2 -> normal
			%yo(:,:,Ch1, Ch2) = norminv(chi2cdf(ynorm, n));
			% just square root for now ...
			yo(:,:,Ch1,Ch2) = sqrt(ynorm);
		else
			% for off-diagonal elements, z transform
			yo(:,:,Ch1, Ch2) = zTrans(y(:,:,Ch1, Ch2));
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