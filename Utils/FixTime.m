function data=FixTime(data)
% Aligns the session to start at time 0 & convert ts from microseconds to
% seconds
%
% spikes
for i=1:numel(data.Spikes)
    data.Spikes{i,1}=data.Spikes{i,1}-data.offset;
%     data.Spikes{i,1}(data.Spikes{i,1}<0)=[];
    % convert microseconds to seconds
    data.Spikes{i,1}=data.Spikes{i,1}./10^6;
end

% lfp
% check to see if already aligned
if ~mean(diff(data.lfp.ts)) == 1/data.lfp.lfpsamplerate || ~data.lfp.ts(data.lfp.lfpsamplerate + 1) == 1
    data.lfp.ts=data.lfp.ts-data.offset;
%     data.lfp.signal(:,data.lfp.ts<0)=[];
%     data.lfp.theta(:,data.lfp.ts<0)=[];
%     data.lfp.theta_phase(:,data.lfp.ts<0)=[];
%     data.lfp.ts(data.lfp.ts<0)=[];
    % convert microseconds to seconds
    data.lfp.ts=data.lfp.ts./10^6;
end

% events
data.events=(data.events-data.offset)./10^6;

% uncorrected frames from linear track
if isfield(data,'linear_track')
    for i=1:length(data.linear_track)
        if isfield(data.linear_track{i},'nonlinearFrames')
            data.linear_track{i}.nonlinearFrames(:,1)=...
                (data.linear_track{i}.nonlinearFrames(:,1)-data.offset)./10^6;
        end
    end
end

%  frames
data.frames(:,1)=(data.frames(:,1)-data.offset)./10^6;

% update structure to reflect change
data.ts_timescale='sec';
end