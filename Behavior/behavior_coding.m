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
p.addParameter('behaviors',{'grooming','rearing','headscan','freezing'});
parse(p,varargin{:});

data.path = p.Results.path;
data.behaviors = p.Results.behaviors;


if isempty(data.path)
    data.path=uigetfile('*.*','find your video you want to analyze');
end

% load video
data.video=VideoReader(data.path);
data.video.CurrentTime=0;

% set up figure
F=figure('Name',[data.path],'NumberTitle','off');
F.Color=[1 1 1];
set(F, 'Position', get(0, 'Screensize'));

subplot('Position',[.05 .05 .7 .90]);
axis off
vidFrame = readFrame(data.video);
image(vidFrame);
axis off

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
    'Callback', @(src,event,data)Done);
playButton = ...
    uicontrol('Parent', F, ...
    'Units', 'Normalized', ...
    'Position', [0.85 0 0.15 0.05], ...
    'Style', 'togglebutton', ...
    'String', 'play', ...
    'Callback', @(src,event)play);
behavior1 = uicontrol('Parent', F, ...
    'Units', 'Normalized', ...
    'Position', [.85 0.95 0.15 0.05], ...
    'Style', 'togglebutton', ...
    'String', data.behaviors{1}, ...
    'Callback', @(src,event)recordtime_behav1);
behavior2 = uicontrol('Parent', F, ...
    'Units', 'Normalized', ...
    'Position', [.85 0.80 0.15 0.05], ...
    'Style', 'togglebutton', ...
    'String', data.behaviors{2}, ...
    'Callback', @(src,event)recordtime_behav2);
behavior3 = uicontrol('Parent', F, ...
    'Units', 'Normalized', ...
    'Position', [.85 0.65 0.15 0.05], ...
    'Style', 'togglebutton', ...
    'String', data.behaviors{3}, ...
    'Callback', @(src,event)recordtime_behav3);
behavior4 = uicontrol('Parent', F, ...
    'Units', 'Normalized', ...
    'Position', [.85 0.50 0.15 0.05], ...
    'Style', 'togglebutton', ...
    'String', data.behaviors{4}, ...
    'Callback', @(src,event)recordtime_behav4);


    function play
        while hasFrame(data.video) && playButton.Value==1
            vidFrame = readFrame(data.video);
            image(vidFrame);
            axis off
            title(['Time ',num2str(data.video.CurrentTime), 'sec'])
            slider_frameToDisplay.Value=data.video.CurrentTime;
            pause(1/data.video.FrameRate);
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
            data.results.(data.behaviors{i}).occurrences=length(data.(data.behaviors{i}))/2;
            
            if ~sum(data.(data.behaviors{i})(:,2))==sum(data.(data.behaviors{i})(:,2)==0)
                error(['Different Number of starts and stops for ',data.behaviors{i}])
            end
            data.results.(data.behaviors{i}).totaltime=sum(diff(data.(data.behaviors{i})(:,1)));
        end
        save_data(data)
        
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
        if ~isfield(data,(data.behaviors{1}))
            data.(data.behaviors{1})(1,:)=[data.video.CurrentTime,behavior1.Value];
            return
        end
        data.(data.behaviors{1})(size(data.(data.behaviors{1}),1)+1,:)=[data.video.CurrentTime,behavior1.Value];
    end

    function recordtime_behav2
        if ~isfield(data,(data.behaviors{2}))
            data.(data.behaviors{2})(1,:)=[data.video.CurrentTime,behavior2.Value];
            return
        end
        data.(data.behaviors{2})(size(data.(data.behaviors{2}),1)+1,:)=[data.video.CurrentTime,behavior2.Value];
    end

    function recordtime_behav3
        if ~isfield(data,(data.behaviors{3}))
            data.(data.behaviors{3})(1,:)=[data.video.CurrentTime,behavior3.Value];
            return
        end
        data.(data.behaviors{3})(size(data.(data.behaviors{3}),1)+1,:)=[data.video.CurrentTime,behavior3.Value];
    end

    function recordtime_behav4
        if ~isfield(data,(data.behaviors{4}))
            data.(data.behaviors{4})(1,:)=[data.video.CurrentTime,behavior4.Value];
            return
        end
        data.(data.behaviors{4})(size(data.(data.behaviors{4}),1)+1,:)=[data.video.CurrentTime,behavior4.Value];
    end

end


