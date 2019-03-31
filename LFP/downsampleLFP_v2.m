function [data]=downsampleLFP_v2(eegfile,data)
warning('off', 'MATLAB:dispatcher:nameConflict')

com=which('downsampleLFP');
com=strsplit(com,filesep);
basedir=[com{1},filesep,'Users',filesep,com{3},filesep,'GoogleDrive',filesep,'MatlabDir'];
addpath(genpath([basedir,filesep,'eeglab14_0_0b']),...
    genpath([basedir,filesep, 'buzcode', filesep, 'externalPackages', filesep, 'FMAToolbox'])...
    ,genpath([basedir,filesep, 'CMBHOME']));

for ii=1:length(eegfile)
    % extract CSC data with mex
    if ismac==1
        [Timestamps,Samples]= Nlx2MatCSC_v3(eegfile{ii}, [1 0 0 0 1], 0, 1);
    else
        [Timestamps,Samples]= Nlx2MatCSC(eegfile{ii}, [1 0 0 0 1], 0, 1);
    end
    
    % reshape Samples variable into a vector
    EEGreshaped = reshape(Samples,1,length(Samples(:)));
        
    if ii==1
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
        
        % preallocate
        EEG_DownSampledData=zeros(length(eegfile),length(EEG_DownSampledTimestamps));
        EEGthetaData=zeros(length(eegfile),length(EEG_DownSampledTimestamps));
        theta_phase=zeros(length(eegfile),length(EEG_DownSampledTimestamps));
        theta_amp=zeros(length(eegfile),length(EEG_DownSampledTimestamps));
        sec=EEG_DownSampledTimestamps/10^6;
        sec=sec-(sec(1));
    end
    
    % check to see if ts and samples has same length (it's rare they won't)
    if length(EEGtsUpSampled)~=length(EEGreshaped)
        [~,I]=max([length(EEGtsUpSampled);length(EEGreshaped)]);
        if I==1
            EEGreshaped=[EEGreshaped,zeros(1,length(EEGtsUpSampled)-length(EEGreshaped))];
        elseif I==2
            EEGreshaped(length(EEGtsUpSampled)+1:end)=[];
        end
    end
    % interpolate data points
    EEG_DownSampledData(ii,:)=interp1(EEGtsUpSampled, EEGreshaped, EEG_DownSampledTimestamps);
    
    %Filter for theta - Takes 'signal' and bandpasses it to theta frequencies (6 to 10 Hz)
    EEGthetaData(ii,:)=ThetaFilter(EEG_DownSampledData(ii,:), NewSFreq);
    %Filter for Delta - Takes 'signal' and bandpasses it for delta frequencies (2 to 4 Hz)
    EEGdeltaData(ii,:) = DeltaFilter(signal, NewSFreq);
    
    %Get theta amplitude
    theta_amp(ii,:) = InstAmplitude(EEGthetaData(ii,:));
    %Get theta phase
    theta_phase(ii,:) = InstPhase(EEGthetaData(ii,:));
    
    
    
    % RUN FMA PHASE CODE
%     [phase,amplitude,~]=Phase([sec',EEGthetaData(ii,:)']);
%     theta_phase(ii,:)=phase(:,2)';
%     theta_amp(ii,:)=amplitude(:,2)';
%     ts_sec=phase(:,1)';
        
    clearvars -except sampleout eegfile i EEG_DownSampledTimestamps EEGtsUpSampled...
        EEG_DownSampledData NewSFreq EEGthetaData sec theta_phase theta_amp ts_sec data
end
% save([extractBefore(eegfile{1,:},'CSC'),'LFP.mat'],'EEG_DownSampledTimestamps',...
%     'EEG_DownSampledData','EEGthetaData','theta_phase','theta_amp','-v7.3')

data.lfp.ts=EEG_DownSampledTimestamps;
data.lfp.signal=EEG_DownSampledData;
data.lfp.theta=EEGthetaData;
data.lfp.theta_phase=theta_phase;
data.lfp.theta_amp=theta_amp;

end