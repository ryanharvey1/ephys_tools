function [yo, fo]=mtcsd1(varargin)
%function [yo, fo]=mtcsd1(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers);
% Multitaper Cross-Spectral Density
% function A=mtcsd1(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers)
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
% series and you'll get a 2-row matrix of cross-spectra
% 
% The output is length( f ) x 2 x nChannels:
%   The first row is the psd of x's columns
%   The second row is the csd the first column and all other columns
%   For instance, yo( :, 2, 3 ) is the csd of x( :, 1 ) and x( :, 3 )
%
% If you want an all-to-all CSD, use mtcsd (n^2 elements; here 2n)
% 
% NB they are cross-spectra not coherences. If you want coherences use
% mtcohere
%
% see also mtcsd, mtcsd2

% Original code by Partha Mitra - modified by Ken Harris
% Also containing elements from specgram.m
% 18-mar-13 ES modified from mtcsd.m

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
nFFTChunks = 1+floor(((nSamples-WinLength)/winstep)); % HOPE THIS WORKS IN GENEREL
% turn this into time, using the sample frequency
t = winstep*(0:(nFFTChunks-1))'/Fs;

% allocate memory now to avoid nasty surprises later
y=complex(zeros(nFFT, 2, nChannels )); % output array
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
    Ch1 = 1;
    Temp1 = squeeze(Periodogram(:,:,Ch1));
    for Ch2 = Ch1:nChannels % don't compute cross-spectra twice
        Temp2 = squeeze(Periodogram(:,:,Ch2));
        Temp2 = conj(Temp2);
        Temp3 = Temp1 .* Temp2;
        eJ=sum(Temp3, 2);
        y(:,2,Ch2)= y(:,2,Ch2) + eJ/(nTapers*nFFTChunks); % cross

        Temp3 = squeeze(Periodogram(:,:,Ch2)) .* Temp2;
        eJ=sum(Temp3, 2);
        y(:,1,Ch2)= y(:,1,Ch2) + eJ/(nTapers*nFFTChunks); % auto
        
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
    for row = 1 : 2
        for Ch2 = 1:nChannels
            subplot(2,nChannels, Ch2+(row-1)*nChannels);
            plot(f,20*log10(abs(y(:,row,Ch2))+eps));
            grid on;
            if(row==2 || Ch2==1)
                ylabel('csd (dB)');
            else
                ylabel('psd (dB)');
            end;
            xlabel('Frequency');
            xlim( f( [ 1 end ] ) )
        end;
    end
elseif nargout == 1
    yo = y;
elseif nargout == 2
    yo = y;
    fo = f;
end

return

% EOF
