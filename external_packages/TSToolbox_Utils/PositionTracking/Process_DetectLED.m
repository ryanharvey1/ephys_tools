% Detect red and blue LEDs position in a video file and creates a 'led' file
%
%  USAGE
%
%    Process_DetectLED(videoFile,<options>)
%
%    videoFile      path to Basler video file, including the '.avi'
%                   extension or not
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'manualROI'   boolean: whether or not you want to manually adjust the
%                   area where LEDs are detected (default = 1)
%    =========================================================================
%
% DEPENDENCIES:
%
%   Computer vision toolbox


% Copyright (C) 2015 Adrien Peyrache, some inspiration from John D Long II
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.


function whl = Process_DetectLED(fbasename,varargin)

%Parameters
manualROI = 0; %Define manually a ROI before LED tracking
smoothVideo = 0; %Apply a Gaussian filter on the image before LED detection
printThres = 0; %Display threshold for LED detection (for debug)


% Parse options
for i = 1:2:length(varargin),
  if ~isa(varargin{i},'char'),
    error(['Parameter ' num2str(i+3) ' is not a property (type ''help Process_ConvertBasler2Pos'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'manualroi',
      manualROI = varargin{i+1};
      if ~isa(manualROI,'numeric')
        error('Incorrect value for property ''manualROI'' (type ''help Process_DetectLED'' for details).');
      end
    case 'smoothvideo',
      smoothVideo = varargin{i+1};
      if ~isa(smoothVideo,'numeric')
        error('Incorrect value for property ''manualROI'' (type ''help Process_DetectLED'' for details).');
      end
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help Process_DetectLED'' for details).']);
  end
end


if strcmp(fbasename(end-3:end),'avi')
    fbasename = fbasename(end-3:end);
end

file = [fbasename '.avi'];

if ~exist(file,'file')
    warning('No video file. Create an empty one')
    whl = [];
    return
end

% Creater readerobj of file
videoObj    = vision.VideoFileReader(file);
videoSize   = videoObj.info.VideoSize;
width       = videoSize(1);
height      = videoSize(2);
threshF     = 0.3;    % minimum threshold on foreground pixel intensity

% Initialize grid for locating centroid of LED
[X,Y] = meshgrid(1:width,1:height);

% Define ROI (manual definition of the environment)
%try 
%First, we adapt the dynamical range of the pixel colors for manual
%selection of the ROI (i.e. the environment)
firstFrame  = step(videoObj);
frame = imadjust(rgb2gray(firstFrame));

if manualROI
    ok = 0;
    figure(1),clf

    while ~ok
        frameNew = frame;
        clf,imshow(frame),caxis([0 0.5])
        fprintf('Define your ROI. Click ''enter'' when finished\n')
        [x,y] = ginput;
        inArea = inpolygon(X(:),Y(:),x,y);
        inArea = double(reshape(inArea,[height width]));
        frameNew(~inArea) = 0;

        clf,imshow(frameNew),caxis([0 0.5])

        reply = input('OK with the result? Y/N [Y]:','s');
        if ~strcmp(reply,'N')
            ok = 1;
        end
    end
    frame = frameNew;
    clear frameNew;
else
    inArea = ones(size(frame));
end


% Initialize background as a grayscale image of the first frame
bg_bw       = rgb2gray(firstFrame);

sqROIx      = any(inArea);
sqROIy      = any(inArea');
height      = sum(sqROIy);
width       = sum(sqROIx);

inArea = reshape(inArea(sqROIy,sqROIx),[height width]);

% Initialize color mask
mask  = zeros(height,width,'uint8');

% Initialize whl matrix
whl = [];
count = 1;

while ~isDone(videoObj)
    
    if mod(count,100)==0
        if count~=100
            backSp = repmat('\b',[1 length(num2str(count-100))]);
            fprintf(backSp)
        end
        fprintf('%i',count)
    end 
    if count~=1
       fr = step(videoObj);
    else
       fr = firstFrame;
    end

    %%% Find centroid of remaining pixels %%%
    %Red against all other
    
    %fr_col        = squeeze(fr(sqROIy,sqROIx,1));% sum(fr(sqROIy,sqROIx,1)));
    fr_col  = fr(sqROIy,sqROIx,3) - max(fr(sqROIy,sqROIx,1:2),[],3);
    fr_col  = squeeze(fr_col);
    fr_col(fr_col<0)  = 0;
    
    if smoothVideo
        fr_col        = gaussFilt(fr_col,1);
    end
    %thresh        = max(percentile(fr_col(fr_col(:)<1),0.999),threshF); %Try to find a good threshold
    thresh = threshF;
    label         = logical(fr_col>thresh & inArea);
    mask(label)   = 1;
    mask(~label)  = 0;
    
    if printThres
        fprintf('thres red: %f\n',thresh)
    end
    
    %bw_mask = squeeze(mask(:,:,1));
    [CC,nc] = bwlabel(mask);
    
    if nc>0
        pixels = regionprops(CC,'PixelList');
        centroidSize = zeros(nc,1);
        for ii=1:nc
            p = pixels(ii).PixelList;
            centroidSize(ii) = length(p);
        end
        [~,mxIx] = max(centroidSize);
        Rr   = round(mean(pixels(mxIx).PixelList,1));
    else
        Rr = [-1 -1];
    end
    
    %Blue against the other colors
   
    %fr_col        = squeeze(fr(sqROIy,sqROIx,3));
    
    fr_col  = fr(sqROIy,sqROIx,1) - max(fr(sqROIy,sqROIx,2:3),[],3);
    fr_col  = squeeze(fr_col);
    fr_col(fr_col<0)  = 0;    
    
    if smoothVideo
        fr_col        = gaussFilt(fr_col,1);
    end
    sortVal    = sort(double(fr_col(fr_col(:)<1)),'ascend');
    thresh        = max(sortVal(round(0.999*length(sortVal))),threshF);
    label         = logical(fr_col>thresh & inArea);
    mask(label)   = 1;
    mask(~label)  = 0;

    if printThres
        fprintf('thres blue: %f\n',thresh)
    end
    
    [CC,nc] = bwlabel(mask);
    pixels = regionprops(CC,'PixelList');
    
    if nc>0
        pixels = regionprops(CC,'PixelList');
        centroidSize = zeros(nc,1);
        for ii=1:nc
            p = pixels(ii).PixelList;
            centroidSize(ii) = length(p);
        end
        [~,mxIx] = max(centroidSize);
        Br   = round(mean(pixels(mxIx).PixelList,1));
    else
        %keyboard
         Br = [-1 -1];
    end
    
     whl = [whl;[Rr(1),Rr(2),Br(1),Br(2)]];
    
    % End processing time
    
    % Display results every 1000 frame 
    if 0 %set it at 1 for debug
        if mod(count,1000)==0
            figure(2),
            plot(whl(:,1),whl(:,2),'r')
            hold on
            plot(whl(:,3),whl(:,4),'b')
            drawnow
            %keyboard
        end
    end
    count = count+1;
end
% catch
%    keyboard
% end

fprintf('Number of undetected points: %d red, %d blue\n\n',sum(whl(:,1)==-1),sum(whl(:,3)==-1))

pos = whl;
pos(pos==-1) = NaN;
if 0
figure(2),clf
    plot(pos(:,1),pos(:,2),'r')
    hold on
    plot(pos(:,3),pos(:,4),'b')
end
%write result file

dlmwrite([fbasename '.led'],whl,'\t');
release(videoObj);