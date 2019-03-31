function TimeFreqPlots(path,StartofRec,EndofRec,event)
%TimeFreqPlots Creates time freq plots
%   Detailed explanation goes here

% extract CSC data with mex
eegfile = FindFiles('*.ncs', 'StartingDirectory',path);

for cscfile=1:length(eegfile)
    if ismac==1
        [Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples]= Nlx2MatCSC_v3(eegfile{cscfile}, [1 1 1 1 1], 0, 1);
    else
        [Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples]= Nlx2MatCSC(eegfile{cscfile}, [1 1 1 1 1], 0, 1);
    end
    % reshape Samples variable into a vector
    EEGreshaped = reshape(Samples,1,length(Samples(:)));
    
    % increase the number of timestamps by a factor of 512 so that there is a timestamps for every data point
    factor = 512;
    EEGtsUpSampled = zeros(factor,length(Timestamps));
    EEGtsUpSampled(1,:)=Timestamps;
    TSdiff=diff(Timestamps);
    TSdiff(end+1)=TSdiff(end);
    TSdiffIncrements=TSdiff/factor;
    for i=1:length(Timestamps)
        for j=2:factor
            EEGtsUpSampled(j,i)=Timestamps(i)+TSdiffIncrements(i)*(j-1);
        end
    end
    EEGtsUpSampled = (EEGtsUpSampled(:))';
    
    % downsample to 1000 Hz (is a little nicer to work with than 1523 Hz (Taube lab) or 2034 Hz (McN lab))
    NewSFreq = 1000;
    
    % calculate number of samples in down sampled data
    DurationData = (EEGtsUpSampled(end)-EEGtsUpSampled(1))/1000000; % in seconds
    NumNSamples = floor(DurationData*NewSFreq);
    
    % create new timestamps
    index = 1:NumNSamples;
    EEG_DownSampledTimestamps(index) = EEGtsUpSampled(1)+(1000000*1/NewSFreq)*(index-1);
    
    % interpolate data points
    EEG_DownSampledData = interp1(EEGtsUpSampled, EEGreshaped, EEG_DownSampledTimestamps);
    
    % RESTRICT DATA BY START AND END OF EVENT
if (length(StartofRec) + length(EndofRec))>2
    EEG_DownSampledData=EEG_DownSampledData(EEG_DownSampledTimestamps>StartofRec(event) & EEG_DownSampledTimestamps<EndofRec(event));
    EEG_DownSampledTimestamps=EEG_DownSampledTimestamps(EEG_DownSampledTimestamps>StartofRec(event) & EEG_DownSampledTimestamps<EndofRec(event));
end

EEG_DownSampledData_ALL(cscfile,:)=EEG_DownSampledData;
end
    % % Compute the multi-taper spectrogram using Chronux toolbox
        figure(1); %subplot(8,1,cscfile);
        dt = 0.001; %Define the time resolution.
        Fs = 1/dt; %The sampling frequency in Hertz
        params.Fs = Fs;
        params.tapers = [2 1]; %Chronux specific parameters
%         num_freqs=50;
%         frex=logspace(.01,1.7,num_freqs);
        [Smtm,Tmtm,Fmtm] = mtspecgramc(zscore(mean(EEG_DownSampledData_ALL)),[1 0.5],params);
         power=zscore(10*log10(Smtm/max(Smtm(:)))');
%         Norm_power=(power-min(min(power)))/max(max(power))-min(min(power));
        h=pcolor(Tmtm, Fmtm, power); axis xy
%         h=pcolor(Tmtm, Fmtm, 10*log10(Smtm/max(Smtm(:)))'); axis xy
        shading interp
%         set(gca,'clim',[-90 0],'ylim',[1 50],'YTick',1:4:length(frex),'YTickLabel',round(frex(1:4:end))); % 'xlim',[0 100]
%         set(gca,'clim',[-90 0],'ylim',[1 50]); % 
        set(gca,'ylim',[1 50]); % 

%          set(gca,'YScale','log')
        colorbar;
        colormap jet
        xlabel('Time (s)')
        ylabel('Freq (Hz)')
end
% end

