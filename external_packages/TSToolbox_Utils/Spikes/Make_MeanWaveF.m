function [spkWidth pk2pk meanWaveF maxIx] = Make_MeanWaveF(fbasename)

% [spkWidth pk2pk meanWaveF maxIx] = Make_MeanWaveF(fbasename)
% 
% this program discriminate interneurons and pyramidal cells 
% It saves two values for each cell: half peak width and peak to peak
% both of them in ms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [pk2Val2 halfPkWidth meanWaveF pk2Val1 peakAsy] = Make_MeanWaveF(fbasename)

xml_data = LoadXml(fbasename);
fs = xml_data.SampleRate;

shankIx = [1:length(xml_data.SpkGrps)];

meanWaveF = {};
spkWidth = [];
pk2pk = [];
maxIx = [];


for ch=shankIx
    ne = length(xml_data.SpkGrps(ch).Channels);
    ns = xml_data.SpkGrps(ch).nSamples;
    
    disp(['Shank #' num2str(ch)])
    fname = [fbasename '.clu.' num2str(ch)];
    if exist(fname,'file')
        clu = load(fname);
        clu = clu(2:end);
        nClu = unique(clu(clu>1));

        for c=1:length(nClu)

            wav = LoadSpikeWaveF([fbasename '.spk.' num2str(ch)], ne,ns,find(clu==nClu(c)));
            w = [];
            for e=1:ne
                w = [w;mean(squeeze(wav(e,:,:)),2)'];
            end
            clear wav
            meanWaveF{end+1} =w ;

            if 0
                figure(1),clf
                for e=1:ne
                    subplot(ne,1,e)
                    plot(w(e,:))
                end
            end

            [dummy,mxIx] = max(sum(w.*w,2));
            maxIx(end+1) = mxIx;

            w = w(mxIx,:)';
            wu = w-mean(w);
            nSamples = length(wu);

            t = [0:1/(fs):(nSamples-1)/(fs)]*1000;

            %We assume that trough position is always 17
            % Trough to peak
            minPos = 17;
            [maxVal2,maxPos2] = max(wu(minPos+1:end));
            mx = maxPos2+minPos;
            p2v = t(mx)-t(minPos);



            [wave f t] = getWavelet(w,20000,900,3000,128);
            %We consider only the central portion of the wavelet because we
            %haven't filtered it before hand (e.g. with a Hanning window)
            wave = wave(:,int16(length(t)/4):3*int16(length(t)/4));

            %Where is the max frequency?
            [maxPow ix] = max(wave);
            [dumy mix] = max(maxPow);
            ix = ix(mix);
            spkW = 1000/f(ix);

            if 0
            t = t(int16(length(t)/4):3*int16(length(t)/4));
            figure(1),clf
            subplot(1,2,1)
                plot(w)
            subplot(1,2,2)
                contourf(t,f,log10(wave))
            disp([p2v spkW])
            pause
            end

           spkWidth(end+1) =spkW;
           pk2pk(end+1) = p2v;
        end
    else
        warning([fname ' doesn''t exist']);
    end
    clear spk;
end