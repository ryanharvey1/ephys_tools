% Clean red and blue LEDs position of a 'led' file
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


% Copyright (C) 2015 Adrien Peyrache, some inspiration from John D Long II
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.


function whl = Process_CleanLED(fbasename,varargin)
manualROI = 1;

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
   
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help Process_DetectLED'' for details).']);
  end
end


if strcmp(fbasename(end-3:end),'led')
    fbasename = fbasename(end-3:end);
end

file = [fbasename '.led'];

if ~exist(file,'file')
    warning('No led file')
    keyboard
end

% Initialize grid for locating centroid of LED

pos = load([fbasename '.led']);
width = max([pos(:,1);pos(:,3)]);
height = max([pos(:,2);pos(:,4)]);

pos(pos == -1) = NaN;

ok = 0;
figure(1),clf

while ~ok
    postmp = pos;
    postmp(postmp == -1) = NaN;

    clf
    plot(pos(:,1),pos(:,2),'r')
    hold on
    plot(pos(:,3),pos(:,4),'b')
    
    fprintf('Define external RED LED ROI. Click ''enter'' when finished\n')
    [xex,yex] = ginput;
    fprintf('Define internal RED LED ROI. Click ''enter'' when finished\n')
    [xin,yin] = ginput;
    inArea = ~inpolygon(pos(:,1),pos(:,2),xex,yex) | inpolygon(pos(:,1),pos(:,2),xin,yin);
    postmp(inArea,1:2) = NaN;
    
    fprintf('Define external BLUE LED ROI. Click ''enter'' when finished\n')
    [xex,yex] = ginput;
    fprintf('Define internal BLUE LED ROI. Click ''enter'' when finished\n')
    [xin,yin] = ginput;
    inArea = ~inpolygon(pos(:,3),pos(:,4),xex,yex) | inpolygon(pos(:,3),pos(:,4),xin,yin);
    postmp(inArea,3:4) = NaN;

    clf
    plot(postmp(:,1),postmp(:,2),'r')
    hold on
    plot(postmp(:,3),postmp(:,4),'b')
    
    reply = input('OK with the result? Y/N [Y]:','s');
    if ~strcmp(reply,'N')
        ok = 1;
    end
end
pos = postmp;
pos(isnan(pos)) = -1;

dlmwrite([fbasename '.led'],pos,'\t')

