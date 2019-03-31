function [yo, fo]=mtcsd(varargin);
%function [yo, fo]=mtcsg(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers);
% Multitaper Cross-Spectral Density
% function A=mtcsd(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers)
% x : input time series
% nFFT = number of points of FFT to calculate (default 1024)
% Fs = sampling frequency (default 2)
% WinLength = length of moving window (default is nFFT)
% nOverlap = overlap between successive windows (default is WinLength/2)
% NW = time bandwidth parameter (e.g. 3 or 4), default 3
% nTapers = number of data tapers kept, default 2*NW -1
%
% output yo is yo(f)
%
% If x is a multicolumn matrix, each column will be treated as a time
% series and you'll get a matrix of cross-spectra out yo(f, Ch1, Ch2)
% NB they are cross-spectra not coherences. If you want coherences use
% mtcohere

% Original code by Partha Mitra - modified by Ken Harris
% Also containing elements from specgram.m

% default arguments and that
[x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers] = mtparam(varargin);
winstep = WinLength - nOverlap;

clear varargin; % since that was taking up most of the memory!

nChannels = size(x, 2);
nSamples = size(x,1);

% check for column vector input
if nSamples == 1
	x = x';
	nSamples = size(x,1);
	nChannels = 1;
end;

% calculate number of FFTChunks per channel
nFFTChunks = round(((nSamples-WinLength)/winstep));
% turn this into time, using the sample frequency
t = winstep*(0:(nFFTChunks-1))'/Fs;

% allocate memory now to avoid nasty surprises later
y=complex(zeros(nFFT, nChannels, nChannels)); % output array
Periodogram = complex(zeros(nFFT, nTapers, nChannels)); % intermediate FFTs
Temp1 = complex(zeros(nFFT, nTapers));
Temp2 = complex(zeros(nFFT, nTapers));
Temp3 = complex(zeros(nFFT, nTapers));
eJ = complex(zeros(nFFT,1));

% calculate Slepian sequences.  Tapers is a matrix of size [WinLength, nTapers]
[Tapers V]=dpss(WinLength,NW,nTapers, 'calc');

% New super duper vectorized alogirthm
% compute tapered periodogram with FFT 
% This involves lots of wrangling with multidimensional arrays.

TaperingArray = repmat(Tapers, [1 1 nChannels]);
for j=1:nFFTChunks
	Segment = x((j-1)*winstep+[1:WinLength], :);
	if (~isempty(Detrend))
		Segment = detrend(Segment, Detrend);
	end;
	SegmentsArray = permute(repmat(Segment, [1 1 nTapers]), [1 3 2]);
	TaperedSegments = TaperingArray .* SegmentsArray;

	Periodogram(:,:,:) = fft(TaperedSegments,nFFT);


	% Now make cross-products of them to fill cross-spectrum matrix
	for Ch1 = 1:nChannels
		for Ch2 = Ch1:nChannels % don't compute cross-spectra twice
			Temp1 = squeeze(Periodogram(:,:,Ch1));
			Temp2 = squeeze(Periodogram(:,:,Ch2));	
			Temp2 = conj(Temp2);
			Temp3 = Temp1 .* Temp2;
			eJ=sum(Temp3, 2);
			y(:,Ch1, Ch2)= y(:,Ch1,Ch2) + eJ/(nTapers*nFFTChunks);

		end
	end
end

% now fill other half of matrix with complex conjugate
for Ch1 = 1:nChannels
	for Ch2 = (Ch1+1):nChannels % don't compute cross-spectra twice
		y(:, Ch2, Ch1) = conj(y(:,Ch1,Ch2));
	end
end
		
% set up f array
if ~any(any(imag(x)))    % x purely real
	if rem(nFFT,2),    % nfft odd
		select = [1:(nFFT+1)/2];
	else
		select = [1:nFFT/2+1];
	end
	y = y(select,:,:,:);
else
	select = 1:nFFT;
end
	
f = (select - 1)'*Fs/nFFT;

% we've now done the computation.  the rest of this code is stolen from
% specgram and just deals with the output stage

if nargout == 0
	% take abs, and plot results
    newplot;
    for Ch1=1:nChannels, for Ch2 = 1:nChannels
    	subplot(nChannels, nChannels, Ch1 + (Ch2-1)*nChannels);
		plot(f,20*log10(abs(y(:,Ch1,Ch2))+eps));
		grid on;
		if(Ch1==Ch2) 
			ylabel('psd (dB)'); 
		else 
			ylabel('csd (dB)'); 
		end;
		xlabel('Frequency');
	end; end;
elseif nargout == 1
    yo = y;
elseif nargout == 2
    yo = y;
    fo = f;
end
