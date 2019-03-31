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
        
        nSamples = length(w);

        t = [0:1/(Fq):(nSamples-1)/(Fq)]*1000;
        %wd = mydetrend(w(:));
        wd = w;
        wd = wd-mean(wd(end-10:end));
        wu = resample(wd,10,1);
        wu = gaussianFilter(wu,10*Fq,[200 6000]);
        Fq = Fq*10;
        
        %Positive or negative spike?
        baseLine = mean(wu(1:7));
        [~,absPk] = max(abs(wu));
        if wu(absPk) > baseLine || absPk>400
            fprintf('Positive spike, skipping\n')
            tr2pk       = NaN;
            spkWidth    = NaN;
            wave        = NaN;
            halfPkW     = NaN;
            figure(1),clf
            plot(wu)
            keyboard
        else

            %Peak to valley, in ms
            [~,minPos] = min(wu);
            try
                [maxVal,maxPos] = max(wu(minPos+1:minPos+200));
            catch
                keyboard
            end
            p2v = 1000*maxPos/Fq;
             
            %Half peak width, in ms
            wu = wu/maxVal;
            th1 = find(wu(minPos:maxPos+minPos)<0.5);
            if any(th1)
                th1 = maxPos-th1(end);
                th2 = find(wu(maxPos+minPos:end)<0.5);
                if any(th2)
                    th2 = th2(1);
                else
                    th2 = th1;
                end
            else
                    th1 = NaN;
                    th2 = NaN;
            end
            hPk = 1000*(th1+th2)/Fq;
            
            [wave,f,t] = getWavelet(wu,Fq,400,5000,128);
            %We consider only the central portion of the wavelet because we
            %haven't filtered it before hand (e.g. with a Hanning window)
%            wave = wave(:,int16(length(t)/4):3*int16(length(t)/4));

            %Where is the max frequency?
           [maxPow ix] = max(wave);
           [dumy mix] = max(maxPow);
           ix = ix(mix);
           spkW = 1000/f(ix);
      

            spkWidth = spkW;
            tr2pk = p2v;
            halfPkW = hPk; 
        end
       
end
