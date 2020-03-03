
classdef detect_epoch
    % detect_epoch: class for the detection of epochs
    
    % pedestal
    % ends_of_track
    % 
    % TO ADD
    % sleep
    
    methods(Static)
        
        function [epoch,epoch_str] = pedestal(data,varargin)
            
            p = inputParser;
            p.addParameter('figs',0); % percent of track length to detect ends
%             p.addParameter('buffer_sec',5); % seconds of time after and before event
            p.parse(varargin{:});
            figs = p.Results.figs;
%             buffer_sec = p.Results.buffer_sec;

%             % get maze boundaries
%             for i = 1:size(data.events,2)
%                 frames = data.frames(data.frames(:,1) > data.events(1,i) & data.frames(:,1) < data.events(2,i),:);
%                 XY = [frames(:,2),frames(:,3)];
%                 XY(any(isnan(XY),2),:) = [];
%                 k = convhull(XY(:,1),XY(:,2));
%                 bound{i} = XY(k,:);
%             end
            
            % locate times between trials
            for i = 1:size(data.events,2) - 1
                events(:,i) = [data.events(2,i);data.events(1,i+1)];
                epoch_str(i) = "inter_trial";
            end
            
            % add before and after session times to events
            if data.events(1,1) ~=0 || data.events(2,end) ~= data.frames(end,1)
                events = [[data.frames(1,1);data.events(1,1)],...
                    events,...
                    [data.events(2,end);data.frames(end,1)]];
                
                epoch_str = ["before_trial",epoch_str,"after_trial"];
                
%                 % duplicate boundaries for each event
%                 bound = [bound;bound];
%                 bound = bound(:);
            end
            
            epoch = events';

            
%             % find when the animal was outside the boundaries of the maze
%             % during pre,post, & intertrial times
%             if figs
%                 figure
%             end
%             for i = 1:size(events,2)
%                 
%                 frames = data.frames(data.frames(:,1) > events(1,i) & data.frames(:,1) < events(2,i),:);
%                 
%                 in = inpolygon(frames(:,2),frames(:,3),bound{i}(:,1),bound{i}(:,2));
%                 
%                 if figs
%                     subplot(1,size(events,2),i)
%                     plot(bound{i}(:,1),bound{i}(:,2),'k')
%                     hold on
%                     plot(frames(:,2),frames(:,3),'.k')
%                     plot(frames(in,2),frames(in,3),'.r')
%                     axis image
%                 end
%                 
%                 first_idx = find(in,1,'first');
%                 last_idx = find(in,1,'last');
%                 
%                 if isempty(last_idx)
%                     epoch(i,:) = events(:,i);
%                 else
%                     epoch(i,:) = [frames(first_idx,1),frames(last_idx,1)];
%                 end
%             end

            if figs
                figure;
                plot(data.frames(:,1),data.frames(:,2))
                ylabel('X')
                xlabel('time')
                title('red indicates inter-trial times')
                hold on
                for i = 1:size(epoch,1)
                    patch([epoch(i,1) epoch(i,1), epoch(i,2) epoch(i,2)],...
                        [min(ylim) max(ylim) max(ylim) min(ylim)], [1 0 0],...
                        'FaceAlpha',.3,'EdgeColor','none')
                end
            end
        end
        
        function epoch = ends_of_track(data,varargin)
            % detect epochs where the animal is at the ends of the track
            
            
            p = inputParser;
            p.addParameter('track_bound',0.15); % percent of track length to detect ends
            p.addParameter('time_bound',1); % time (sec) where the animal is in boundary
            p.parse(varargin{:});
            track_bound = p.Results.track_bound;
            time_bound = p.Results.time_bound;
            
            track_idx = contains(data.mazetypes,'track','IgnoreCase',true);
            
            events = data.events(:,track_idx);
            
            for i = 1:size(events,2)
                frames = data.frames(data.frames(:,1) > events(1,i) & data.frames(:,1) < events(2,i),:);
                
                bound = range(frames(:,2)) * track_bound;
                
                idx = frames(:,2) < min(frames(:,2)) + bound |...
                    frames(:,2) > max(frames(:,2)) - bound;
                
                % in range for at least 1 second
                [idx]=contiguousframes(idx,data.samplerate * time_bound);
                
                % locate groups
                [start,ends,~]=findgroups(idx);
                
                epoch{i} = [start',ends'];
                
                % figure for testing
                %                 current_epoch = frames(idx,:);
                %                 figure;
                %                 raw_frames = data.linear_track{1, 1}.nonlinearFrames;
                %                 plot(raw_frames(:,2),raw_frames(:,3),'.k')
                %                 hold on
                %                 x = interp1(raw_frames(:,1), raw_frames(:,2), current_epoch(:,1));
                %                 y = interp1(raw_frames(:,1), raw_frames(:,3), current_epoch(:,1));
                %                 plot(x,y,'.r')
                %                 axis image
                
                
            end
        end
    end
end