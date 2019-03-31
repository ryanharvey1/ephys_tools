function [yo, fo, to]=mtpsg(varargin);
%function [yo, fo, to]=mtcsg(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers);
% Multitaper Time-Frequency Power-Spectrum (power spectrogram)
% function A=mtpsg(x,nFFT,Fs,WinLength,nOverlap,NW,nTapers)
% x : input time series
% nFFT = number of points of FFT to calculate (default 1024)
% Fs = sampling frequency (default 2)
% WinLength = length of moving window (default is nFFT)
% nOverlap = overlap between successive windows (default is WinLength/2)
% NW = time bandwidth parameter (e.g. 3 or 4), default 3
% nTapers = number of data tapers kept, default 2*NW -1
%
% output yo is yo(f, t)
%
% The difference between this and mtcsg is that this one computes power
% spectra only, no cross-spectra.  So If x is a multicolumn matrix, 
% each column will be treated as a time series and you'll get a 3D matrix 
% of power-spectra out yo(f, t, Ch)


[x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,nChannels,nSamples,nFFTChunks,winstep,select,nFreqBins,f,t] = mtparam(varargin);

% allocate memory now to avoid nasty surprises later
y=complex(zeros(nFreqBins, nFFTChunks, nChannels)); % output array

for Ch=1:nChannels
	[y(:,:,Ch) f t] = mtcsg(x(:,Ch),nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers);
end

% we've now done the computation.  the rest of this code is stolen from
% specgram and just deals with the output stage

if nargout == 0
	% take abs, and use image to display results
    newplot;
    for Ch=1:nChannels
    	subplot(1, nChannels, Ch);
	    if length(t)==1
	        imagesc([0 1/f(2)],f,20*log10(abs(y(:,:,Ch))+eps));axis xy; colormap(jet)
	    else
	        imagesc(t,f,20*log10(abs(y(:,:,Ch))+eps));axis xy; colormap(jet)
	    end
	end;
    xlabel('Time')
    ylabel('Frequency')
elseif nargout == 1,
    yo = y;
elseif nargout == 2,
    yo = y;
    fo = f;
elseif nargout == 3,
    yo = y;
    fo = f;
    to = t;
end

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu