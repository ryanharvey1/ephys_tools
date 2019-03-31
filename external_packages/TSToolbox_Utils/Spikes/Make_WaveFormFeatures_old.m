function [spkWidth,tr2pk,halfPkW,maxIx,wave] = Make_WaveFormFeatures(waveforms,varargin)


%  Extract spike width, trough to peak and peak half-width from waveform(s)
%  The program assumes spikes are negative. Exclude positive spikes for
%  now.
%
%  USAGE
%
%    Make_WaveFormFeatures(waveforms,<Fq>)
%
%    INPUT
%    waveforms      a vector or matrix of average waveforms (if recorded on a
%                   polytrode). Matrix should be Channels X Samples
%                   OR a cell array of matrices for different neurons
%    Fq (optional)  sampling frequency (default 20,000Hz)
%
%    OUTPUT
%    spkWidth       total spike width (inverse of peak frequency in a
%                   wavelet transform)
%    tr2pk          trough to peak
%    maxIx          channel index of max spike amplitude

% Copyright (C) 2017 Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

Fq = 20000;

if ~isempty(varargin)
   Fq = varargin{1};
end


spkWidth = [];
tr2pk = [];
maxIx = [];
wave = [];
halfPkW  = [];

if isa(waveforms,'cell')%call itself for each item of the cell array
    
    for c=1:length(waveforms)
        [spkW,t2p,halfW,mxI,wv]   = Make_WaveFormFeatures(waveforms{c},Fq);
        spkWidth            = [spkWidth;spkW];
        tr2pk               = [tr2pk;t2p];
        maxIx               = [maxIx;mxI];
        halfPkW             = [halfPkW;halfW];
        wave(end+1,:,:)     = wv;
    end
    
else
        %if more than one channel, select the one with highest amplitude
        %Restrict to the 1sst 32 samples as sometimes late samples
        %contribute to the power
        if min(size(waveforms))>1
            l               = min(32,size(waveforms,2));
            [~,mxIx]        = max(sum(waveforms(:,1:l).^2,2));
            maxIx(end+1)    = mxIx;
            w               = waveforms(mxIx,:);
        else
            w               = waveforms;
        end
        
        wu = w-mean(w);
        nSamples = length(wu);

        t = [0:1/(Fq):(nSamples-1)/(Fq)]*1000;

        %Positive or negative spike?
        baseLine = mean(wu(1:7));
        [~,absPk] = max(abs(wu));
        if wu(absPk) > baseLine || absPk>40
            fprintf('Positive spike, skipping\n')
            tr2pk       = NaN;
            spkWidth    = NaN;
            wave        = NaN;
            halfPkW     = NaN;
        else
        
            % Trough is supposed to be sample #17, but we never know...
            [~,minPos] = min(wu); %and should be the same as absPk...

            % Where is the following peak? Restrict to the 1st 24 samples
            l = min(24,size(waveforms,2)-minPos);
            [maxVal,maxPos] = max(wu(minPos+1:minPos+l));
            maxPos = maxPos+minPos;
            p2v = t(maxPos)-t(minPos);
            
            %[wave,f] = cwt(wu,Fq,'VoicesPerOctave',48);
            %wave = wave(f>500 & f<3000,:);
            %f = f(f>500 & f<3000);

            %Where is the max power?
            %[maxPow,maxF] = max(wave);
            %Which frequency does it correspond to?
            %[~,fIx] = max(maxPow);
            
            %maxF = maxF(fIx);

            %Spike width is the inverse of the peak frequency
            %spkW = 1000/f(maxF);
            spkW = NaN;
            
            %Half Peak Width
            %wu = wu(:).*hanning(length(wu));
            baseLine = mean(wu(end-5:end));
            [maxVal,maxPos] = max(wu(minPos+1:minPos+l));
            maxPos = maxPos+minPos-1;
            %try
            ix1 = minPos;
            if maxVal<0
                halfW = NaN;
            else
                
                while wu(ix1)<maxVal/2 &ix1<length(wu)
                    ix1=ix1+1;
                end
                ix1=ix1-1;
                ix2 = LocalMinima(wu(ix1:end),5,baseLine);
                
                if isempty(ix2)

                    ix2 = ix1+1;
                    while wu(ix2)>maxVal/2 & ix2<length(wu)
                        ix2=ix2+1;
                    end
                else
                    ix2 = ix1+ix2(1)-1;
                    halfVal = (wu(maxPos) + wu(ix2))/2;
                    ix2 = ix1;
                    while wu(ix2)>halfVal & ix2<length(wu)
                        ix2=ix2+1;
                    end
                end
                
%                 if ix2 == length(wu)
%                     halfW = NaN;
%                 else
                    halfW = t(ix2)-t(ix1);
%                 end
            end
            %catch
            %    
             %   figure(1),clf
             %   plot(wu)
             %   keyboard
            %end
            if p2v<0.4 &  halfW>0.3
                figure(1),clf
                plot(w)
                keyboard
            end
            spkWidth =spkW;
            tr2pk = p2v;
            halfPkW = halfW;
        end
       
end
