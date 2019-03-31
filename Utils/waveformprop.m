function prop=waveformprop(waves)
% waveformprop: calculates features on average waveforms
%
% Input:
%           waves: nested cell array number of cells long which includeds 
%           4 wire x 150us matrix for each cell
%
% Output:
%           prop: data structure containing the following features
%                   spikewidth: just time of valley - time of peak. 
%                   ratio of peak height:valley depth
%
% Ryan E Harvey 2018

% cycle through nested cell array

% factor to convert bined time to ms
To_ms=1/length(waves{1});

for i=1:length(waves)
    
    % current wave
    curwave=waves{i};
    
    % find largest wave (wire closest to the cell)
    [p,I]=max(max(curwave,[],2));
    % extract wave 
    peakwave=curwave(I,:);
  
    %% halfwidth
    [p,peakloc]=max(peakwave);
    [~,locs,widths,~]=findpeaks(peakwave);
    halfwidth=widths(locs==peakloc)*To_ms;
    if isempty(halfwidth)
        halfwidth=NaN;
    end
    
    
    %% peak
    [peak, ipeak] = max(max(curwave));
    % valley after peak
    tempcurwave=[zeros(size(curwave,1),ipeak-1)+peak,curwave(:,[ipeak:end])];
    [vlly, ivlly] = min(min(tempcurwave));

    
    %% spike width
    spikewidth=(ivlly-ipeak)*To_ms;
    
    %% peak to valley ratio
    peak2valleyratio=abs(peak)./abs(vlly);
    
    %% Slope
    peak2valleyslope=(vlly-peak)./(ivlly-ipeak);
    
    %% Spike Amplitude
    SpikeAmplitude=abs(peak)-abs(vlly);
    
    %% collect properties
    prop(i,:)=[peak,vlly,spikewidth,halfwidth,peak2valleyratio,peak2valleyslope,SpikeAmplitude];
end

%     w = MaxWaves(:,a);
%     w = [w(1)*ones(1000,1);w;w(end)*ones(1000,1)];
% waves(1,:)
%     [wave f t] = getWavelet(waves(1,:)',32000,500,3000,150);

%       [wave,f,t]=getWavelet(maxmeanwv',20000,500,3000,128,0,0);
%   meanWavelet=mean(wave,1);
%    figure;
%    imagesc(t,f,wave);
%    troughWavelet=wave(:,troughPos);
%    set(gca,'YDir','normal');
%   SpikeFreq=f(meanWavelet==max(meanWavelet));
%   invF=(1/SpikeFreq)*1000;%(= spike width in ms)

end

