function  behavior_coding(varargin)
% behavior_coding: read in video and label different behaviors
%
% Input:    
%           path: path to video file (uigetfile by defalt)
%           behaviors: cell array of 4 behaviors to label (see default below)
%
% output: saves a .mat file containing a structure of occurrences of behaviors 
%           to the same path as video
%
%
% how to use: you can play and pause the video with "play" as well as
% scroll through the video with the scroll bar. When you see a behavior,
% toggle the respective behavior button on (it will change color). When the
% behavior stops, toggle the button off. When done, hit done and the
% results file will be saved to your directory. 
%
% Future work: 
%       make everything robust to user error. 
%       add an option for more behaviors
%       check if multiple video formats work
%       clean up UI
%       make button color change more obvious

% ryan harvey 2019

% parse input
p = inputParser;
p.addParameter('path',[]);
p.addParameter('behaviors',{'grooming','rearing','head_scan','Trial_Start','Trial_End'});
parse(p,varargin{:});

clear global
global data

data.path = p.Results.path;
data.behaviors = p.Results.behaviors;


if isempty(data.path)
    data.path=uigetfile('*.*','find your video you want to analyze');
end

% load video
data.video=VideoReader(data.path);
data.video.CurrentTime=0;
data.playbackspeed=1;

% set up figure
F=figure('Name','subject','NumberTitle','off');
F.Color=[1 1 1];
set(F, 'Position', get(0, 'Screensize'));

subplot('Position',[.05 .05 .7 .90]);
axis off
vidFrame = readFrame(data.video);
image(vidFrame);
axis off
darkBackground(F,[0.2 0.2 0.2],[0.9 0.7 0.7])


slider_frameToDisplay = ...
    uicontrol('Parent', F, ...
    'Units','Normalized', ...
    'Position', [0.05 0 0.7 0.05], ...
    'Style', 'slider', ...
    'SliderStep', [1/data.video.duration, 10/data.video.duration], ...
    'min', 0, 'max', data.video.duration, ...
    'value', 0, ...
    'Callback', @(src,event)slide_vid);
DoneButton = ...
    uicontrol('Parent', F, ...
    'Units', 'Normalized', ...
    'Position', [0 .95 0.15 0.05], ...
    'Style', 'PushButton', ...
    'String', 'DONE (save)', ...
    'FontSize',15,...
    'FontWeight','bold',...
    'Callback', @(src,event,data)Done);
playButton = ...
    uicontrol('Parent', F, ...
    'Units', 'Normalized', ...
    'Position', [0.85 0 0.15 0.05], ...
    'Style', 'togglebutton', ...
    'String', 'start/stop', ...
    'FontSize',15,...
    'FontWeight','bold',...
    'Callback', @(src,event)play);
playfastButton = ...
    uicontrol('Parent', F, ...
    'Units', 'Normalized', ...
    'Position', [0.85 .15 0.15 0.05], ...
    'Style', 'togglebutton', ...
    'String', 'fast play', ...
    'FontSize',15,...
    'FontWeight','bold',...
    'Callback', @(src,event)playfast);
behavior1 = uicontrol('Parent', F, ...
    'Units', 'Normalized', ...
    'Position', [.85 0.95 0.15 0.05], ...
    'Style', 'togglebutton', ...
    'String', data.behaviors{1}, ...
    'FontSize',15,...
    'FontWeight','bold',...
    'Callback', @(src,event)recordtime_behav1);
behavior2 = uicontrol('Parent', F, ...
    'Units', 'Normalized', ...
    'Position', [.85 0.80 0.15 0.05], ...
    'Style', 'togglebutton', ...
    'String', data.behaviors{2}, ...
    'FontSize',15,...
    'FontWeight','bold',...
    'Callback', @(src,event)recordtime_behav2);
behavior3 = uicontrol('Parent', F, ...
    'Units', 'Normalized', ...
    'Position', [.85 0.65 0.15 0.05], ...
    'Style', 'togglebutton', ...
    'String', data.behaviors{3}, ...
    'FontSize',15,...
    'FontWeight','bold',...
    'Callback', @(src,event)recordtime_behav3);
behavior4 = uicontrol('Parent', F, ...
    'Units', 'Normalized', ...
    'Position', [.85 0.50 0.15 0.05], ...
    'Style', 'togglebutton', ...
    'String', data.behaviors{4}, ...
    'FontSize',15,...
    'FontWeight','bold',...
    'Callback', @(src,event)recordtime_behav4);
behavior5 = uicontrol('Parent', F, ...
    'Units', 'Normalized', ...
    'Position', [.85 0.35 0.15 0.05], ...
    'Style', 'togglebutton', ...
    'String', data.behaviors{5}, ...
    'FontSize',15,...
    'FontWeight','bold',...
    'Callback', @(src,event)recordtime_behav5);

    function play
        if playButton.Value==1
            playButton.BackgroundColor='r';
        else
            playButton.BackgroundColor=[0.9400 0.9400 0.9400];
        end
        while hasFrame(data.video) && playButton.Value==1
            vidFrame = readFrame(data.video);
            image(vidFrame);
            axis off
            title(['Time ',num2str(data.video.CurrentTime), 'sec'])
            slider_frameToDisplay.Value=data.video.CurrentTime;
            pause((1/data.video.FrameRate)/data.playbackspeed);
        end
    end

    function playfast
        if playfastButton.Value==1
            data.playbackspeed=5; %was 5, changed by LB 09/june/2019. 5 is max without causing delay. 
            playfastButton.BackgroundColor='r';
        else
            data.playbackspeed=3; %was 1, changed by LB 09/june/2019
            playfastButton.BackgroundColor=[0.9400 0.9400 0.9400];
            
        end
    end

    function slide_vid
        data.video.CurrentTime=slider_frameToDisplay.Value;
        vidFrame = readFrame(data.video);
        image(vidFrame);
        axis off
        title(['Time ',num2str(data.video.CurrentTime), 'sec'])
    end

    function Done
        % check and store data
        for i=1:length(data.behaviors)
            if isfield(data,(data.behaviors{i}))
                data.results.(data.behaviors{i}).occurrences=length(data.(data.behaviors{i}))/2;
                
                if ~sum(data.(data.behaviors{i})(:,2))==sum(data.(data.behaviors{i})(:,2)==0)
                    error(['Different Number of starts and stops for ',data.behaviors{i}])
                end
                data.results.(data.behaviors{i}).totaltime=sum(diff(data.(data.behaviors{i})(:,1)));
            end
        end
        save_data(data)
        
        clear data
        clear global
        disp('closing main fig')
        close(F);
    end

    function save_data(data)
        processedpath=strsplit(data.path,filesep);
        filename=processedpath(end);
        processedpath(end)=[];
        savepath=fullfile(strjoin(processedpath,filesep),extractBefore(filename,'.mov'));
        save(savepath{1},'-struct','data','-v7.3')
        
        disp(['saving data to ',savepath{1}])
    end

    function recordtime_behav1
        if behavior1.Value==1
            behavior1.BackgroundColor='r';
        else
            behavior1.BackgroundColor=[0.9400 0.9400 0.9400];
        end
        
        if ~isfield(data,(data.behaviors{1}))
            data.(data.behaviors{1})(1,:)=[data.video.CurrentTime,behavior1.Value];
            return
        end
        data.(data.behaviors{1})(size(data.(data.behaviors{1}),1)+1,:)=[data.video.CurrentTime,behavior1.Value];
    end

    function recordtime_behav2
        if behavior2.Value==1
            behavior2.BackgroundColor='r';
        else
            behavior2.BackgroundColor=[0.9400 0.9400 0.9400];
        end
        if ~isfield(data,(data.behaviors{2}))
            data.(data.behaviors{2})(1,:)=[data.video.CurrentTime,behavior2.Value];
            return
        end
        data.(data.behaviors{2})(size(data.(data.behaviors{2}),1)+1,:)=[data.video.CurrentTime,behavior2.Value];
    end

    function recordtime_behav3
        if behavior3.Value==1
            behavior3.BackgroundColor='r';
        else
            behavior3.BackgroundColor=[0.9400 0.9400 0.9400];
        end
        if ~isfield(data,(data.behaviors{3}))
            data.(data.behaviors{3})(1,:)=[data.video.CurrentTime,behavior3.Value];
            return
        end
        data.(data.behaviors{3})(size(data.(data.behaviors{3}),1)+1,:)=[data.video.CurrentTime,behavior3.Value];
    end

    function recordtime_behav4
        if behavior4.Value==1
            behavior4.BackgroundColor='r';
        else
            behavior4.BackgroundColor=[0.9400 0.9400 0.9400];
        end
        if ~isfield(data,(data.behaviors{4}))
            data.(data.behaviors{4})(1,:)=[data.video.CurrentTime,behavior4.Value];
            return
        end
        data.(data.behaviors{4})(size(data.(data.behaviors{4}),1)+1,:)=[data.video.CurrentTime,behavior4.Value];
    end

    function recordtime_behav5
        if behavior5.Value==1
            behavior5.BackgroundColor='r';
        else
            behavior5.BackgroundColor=[0.9400 0.9400 0.9400];
        end
        if ~isfield(data,(data.behaviors{5}))
            data.(data.behaviors{5})(1,:)=[data.video.CurrentTime,behavior5.Value];
            return
        end
        data.(data.behaviors{5})(size(data.(data.behaviors{5}),1)+1,:)=[data.video.CurrentTime,behavior5.Value];
    end
end


