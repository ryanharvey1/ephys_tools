function [PDD_max,PDD_min, meanHD,N,plotIdx] = HDdrift(ExtractedAngle,data_video_spk,data_video_nospk,fig)
%Calculates the amount of preferred direction drift via Yoder & Taube 2009.
%By: LB & RH 2018
%
%   INPUT:
%       HeadDirection: Vector of instantaneous head direction
%       FiringRate:    Vector of instantaneous firing rate
%  OUTPUT:
%       PDD: preferred direction drift
% drift_cat: drift category (<360% drift or >360% drift)
if sum(data_video_spk(:,6))==0 %Check for spikes if none, outputs are NAN
    PDD_max=NaN;
    PDD_min=NaN;
    meanHD=NaN;
    N=NaN;
    plotIdx=NaN;
else
    
    bin=1;
    for i=1:5:length(ExtractedAngle)
        if i>=length(ExtractedAngle)-4
            meanHD(bin,1)=circ_mean(ExtractedAngle(i:end,:));
        else
            meanHD(bin,1)=circ_mean(ExtractedAngle(i:i+4,:));
        end
        bin=bin+1;
    end
    clear bin
    
    % bin over 5 frames
    spks=data_video_spk(data_video_spk(:,6)==1,1);
    edges=data_video_nospk(1:5:end,1);
    [N,~] = histcounts(spks,edges);
    % convert to spikes/sec
    N=N*(30/5);
    % smooth over 200ms
    N=smooth(N,6);
    
    plotIdx=N>=max(N)*.5; %find bins where the firing rate is %50 of the peak bin
    if fig==1
        figure; plot(meanHD,'-','Color',[.75 .75 .75],'LineWidth',3); hold on;
        scatter(1:length(meanHD),meanHD(:,:),25,[.5 .5 .5],'o','filled'); hold on;
        y=1:length(meanHD);
        scatter(y(plotIdx==1),meanHD(plotIdx==1,:),50,'r','o','filled');
    elseif fig==0
    end
    
    PDD_min= min(abs(diff(meanHD(plotIdx==1)))); %range of the absolute deviation of PD shifts across entire session
    PDD_max= max(abs(diff(meanHD(plotIdx==1)))); %range of the absolute deviation of PD shifts across entire session

    if isempty(PDD_min) || isempty(PDD_max)
        PDD_min=NaN;
        PDD_max=NaN;
    end
end
end

