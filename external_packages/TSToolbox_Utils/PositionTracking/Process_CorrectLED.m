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

ok = 0;
figure(1),clf

pos(pos == -1) = NaN;

while ~ok

    clf
    plot(pos(:,1),'r')
    hold on
    plot(pos(:,3),'b')
    
    fprintf('Define X range to correct. Click ''enter'' when finished\n')
    [x,y] = ginput;
    
    offset = input('How much offset? Y/N [Y]:','s');
    offset = str2num(offset);
    x(1) = max(x(1),1);
    x(2) = min(x(2),size(pos,1));
    
    ix = round(x(1)):round(x(2));
    pos(ix,1) = pos(ix,1)+offset;
    pos(ix,3) = pos(ix,3)+offset;
    
    reply = input('Continue? Y/N [Y]:','s');
    if strcmp(reply,'N')
        ok = 1;
    end
end

ok = 0;

while ~ok

    clf
    plot(pos(:,2),'r')
    hold on
    plot(pos(:,4),'b')
    
    fprintf('Define Y range to correct. Click ''enter'' when finished\n')
    [x,y] = ginput;
    
    offset = input('How much offset? Y/N [Y]:','s');
    offset = str2num(offset);
    ix = round(x(1)):round(x(2));
    x(1) = max(x(1),1);
    x(2) = min(x(2),size(pos,1));
    pos(ix,2) = pos(ix,2)+offset;
    pos(ix,4) = pos(ix,4)+offset;
    
    reply = input('Continue? Y/N [Y]:','s');
    if strcmp(reply,'N')
        ok = 1;
    end
end

pos(isnan(pos)) = -1;
dlmwrite([fbasename '.led'],pos,'\t')

