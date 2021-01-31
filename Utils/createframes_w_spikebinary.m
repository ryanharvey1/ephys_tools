function [data_video_spk,data_video_nospk]=createframes_w_spikebinary(data,event,cell)
% createframes_w_spikebinary: interpolates position data from spike times, velocity filters
% out positions and spikes (<3cm/second for 1 second)
%
% 
data_video=data.frames;
data_video=data_video(data_video(:,1)>data.events(1,event)...
    & data_video(:,1)<data.events(2,event),:);

SpikeFile=data.Spikes{cell};

% RESTRICT DATA BY START AND END OF EVENT
SpikeFile=SpikeFile(SpikeFile(:,1)>data.events(1,event)...
    & SpikeFile(:,1)<data.events(2,event),:);

% FIND INDEX FOR VELOCITY FILTER (<3cm/second for 1 frame)
in=contiguousframes(data_video(:,5)<3,1);
data_video_nospk=[data_video,in];

% INTERPOLATE SPIKES TO TIMESTAMPS, POSITION, AND VEL DATA
[~,idx] = unique(data_video_nospk(:,1));

X=interp1(data_video_nospk(idx,1),data_video_nospk(idx,2),SpikeFile,'linear');
Y=interp1(data_video_nospk(idx,1),data_video_nospk(idx,3),SpikeFile,'linear');
A = circular_interp(data_video_nospk(idx,1),data_video_nospk(idx,4),SpikeFile);
% A=interp1(data_video_nospk(idx,1),data_video_nospk(idx,4),SpikeFile,'nearest');
VEL=interp1(data_video_nospk(idx,1),data_video_nospk(idx,5),SpikeFile,'linear');
VELidx=interp1(data_video_nospk(idx,1),data_video_nospk(idx,6),SpikeFile,'linear');

% CONCAT AND SORT
data_video_spk=sortrows([[SpikeFile X Y A VEL VELidx ones(size(SpikeFile,1),1)];...
    [data_video_nospk,zeros(length(data_video_nospk),1)]],1);
% clear TS X Y A VEL VELidx

% VELO FILTER BASED ON INDEX CREATED ABOVE
data_video_nospk(logical(in),:)=[];
data_video_nospk(:,6)=[];
data_video_spk(data_video_spk(:,6)==1,:)=[];
data_video_spk(:,6)=[];
end