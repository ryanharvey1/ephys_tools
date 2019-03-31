% ConvertBasler2Whl_All - Convert Basler data of all the sessions and
% concatenate them
%
%  USAGE
%
%    ConvertBasler2Whl_All(filebasename,<options>)
%
%    filebasename   basename of the video file (should be filebasenmae.avi)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%    =========================================================================
%
% DEPENDENCIES:
%
%   LoadBinary
%   Process_DetectLED
%   Process_ConvertBasler2Pos


% Copyright (C) 2015 Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.


function Process_ConvertBasler2Pos_All(fbasename)

% try
    
    folders = dir([fbasename '-*']);
    nRec = length(folders);

    recName = {};
    startTime = [];
    for ii=1:nRec

        if folders(ii).isdir;
            fname = folders(ii).name;
            recName = [recName;{fname}];
            k = strfind(fname,'-');
            k = k(end);
            startTime = [startTime;str2num(fname(k+1:end))];
        end

    end
    [startTime,ix] = sort(startTime);
    recName = recName(ix);
    nRec = length(recName);
    
    posAll = [];
    recDuration = 0;
    
    for ii=1:nRec
        disp(recName{ii})
        cd(recName{ii});
        
        if 1% ~exist([recName{ii} '.pos'],'file')
            
            Process_ConvertBasler2Pos(recName{ii});
        end
        pos = load([recName{ii} '.pos']);
        posAll = [posAll;[pos(:,1)+recDuration pos(:,2:end)]];
        
        if ii~=nRec
            dat = LoadBinary([recName{ii} '_digitalin.dat']);
            recDuration = recDuration+length(dat)/20000;
        end
        cd('..')
    end
        
    fbasename = fullfile(fbasename,fbasename);
    dlmwrite([fbasename '.pos'],posAll,'delimiter','\t','precision',9);
    
% catch
%    keyboard
% end
end

