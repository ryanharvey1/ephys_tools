function [yo, fo]=mtcsd2(x1,varargin)
%function [yo, fo]=mtcsd2(x,y,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers);
% Multitaper Cross-Spectral Density
% function A=mtcsd2(x,y,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers)
% x,y : input time series
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
% series and you'll get a 3-row matrix of cross-spectra
%
% The output is length( f ) x 3 x nColumns:
%   The first row is the psd of x's columns
%   The second row is the psd of y's columns
%   The third row is the csd of x-y
%   For instance, yo( :, 3, 2 ) is the csd of x( :, 2 ) and y( :, 2 )
% 
% If you want an all-to-all CSD, use mtcsd (4n^2 elements; here 6n)
%
% NB they are cross-spectra not coherences. If you want coherences use
% mtcohere
%
% see also mtcsd, mtcsd2

% Original code by Partha Mitra - modified by Ken Harris
% Also containing elements from specgram.m
% 18-mar-13 ES modified from mtcsd.m

% default arguments and that
[x2,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers] = mtparam(varargin);
winstep = WinLength - nOverlap;
if ~isequal( size( x1 ), size( x2 ) )
    error( 'mtcsd2 requires two identical 2D arrays at the input' )
end
clear varargin; % since that was taking up most of the memory!

nColumns = size(x1, 2); % actually, ntrials (there are 2 channels, x1 and x2)
nSamples = size(x1,1);

% check for column vector input
if nSamples == 1
	x1 = x1';
    x2 = x2';
	nSamples = size(x1,1);
	nColumns = 1;
end;

% calculate number of FFTChunks per channel
nFFTChunks = 1+floor(((nSamples-WinLength)/winstep)); % HOPE THIS WORKS IN GENEREL
% turn this into time, using the sample frequency
t = winstep*(0:(nFFTChunks-1))'/Fs;

% allocate memory now to avoid nasty surprises later
y=complex(zeros(nFFT, 3, nColumns)); % output array
Periodogram = complex(zeros(nFFT, nTapers, nColumns)); % intermediate FFTs
Temp1 = complex(zeros(nFFT, nTapers));
Temp2 = complex(zeros(nFFT, nTapers));
Temp3 = complex(zeros(nFFT, nTapers));
eJ = complex(zeros(nFFT,1));

% calculate Slepian sequences.  Tapers is a matrix of size [WinLength, nTapers]
[Tapers V]=dpss(WinLength,NW,nTapers, 'calc');

% New super duper vectorized alogirthm
% compute tapered periodogram with FFT 
% This involves lots of wrangling with multidimensional arrays.

TaperingArray = repmat(Tapers, [1 1 nColumns]);
for j=1:nFFTChunks
	Segment1 = x1((j-1)*winstep+[1:WinLength], :);
	Segment2 = x2((j-1)*winstep+[1:WinLength], :);
	if (~isempty(Detrend))
		Segment1 = detrend(Segment1, Detrend);
		Segment2 = detrend(Segment2, Detrend);
	end;
	SegmentsArray1 = permute(repmat(Segment1, [1 1 nTapers]), [1 3 2]);
	SegmentsArray2 = permute(repmat(Segment2, [1 1 nTapers]), [1 3 2]);
	TaperedSegments1 = TaperingArray .* SegmentsArray1;
	TaperedSegments2 = TaperingArray .* SegmentsArray2;

	Periodogram1(:,:,:) = fft(TaperedSegments1,nFFT);
	Periodogram2(:,:,:) = fft(TaperedSegments2,nFFT);

	% Now make cross-products of them to fill cross-spectrum matrix
    for Col = 1:nColumns
        Temp1 = squeeze(Periodogram1(:,:,Col));
        Temp2 = squeeze(Periodogram2(:,:,Col));
        Temp2 = conj(Temp2);
        Temp3 = Temp1 .* Temp2;
        eJ=sum(Temp3, 2);
        y(:,3, Col)= y(:,3,Col) + eJ/(nTapers*nFFTChunks); % cross
        
        Temp2 = conj(Temp1);
        Temp3 = Temp1 .* Temp2;
        eJ=sum(Temp3, 2);
        y(:,1, Col)= y(:,1,Col) + eJ/(nTapers*nFFTChunks); % auto1
        
        Temp1 = squeeze(Periodogram2(:,:,Col));
        Temp2 = conj(Temp1);
        Temp3 = Temp1 .* Temp2;
        eJ=sum(Temp3, 2);
        y(:,2, Col)= y(:,2,Col) + eJ/(nTapers*nFFTChunks); % auto2
    end
end

% set up f array
if ~any(any(imag(x1))) && ~any(any(imag(x2)))    % x purely real
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
    for row = 1 : 3
        for Col = 1:nColumns
            subplot(2,nColumns, Col+(row-1)*nColumns);
            plot(f,20*log10(abs(y(:,row,Col))+eps));
            grid on;
            if(row==3 || Col==1)
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

