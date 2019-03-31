% ConvertBasler2Whl - Convert Basler data to whl file format
%
%  USAGE
%
%    ConvertBasler2Whl(videoFile,digitalSig,<options>)
%
%    videoFile      path to Basler video file
%    digitalSig     path to dat file containing the syn TTLs
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'nChannels'   number of data channels in the file (default = 1)
%     'channels'    channels to read (default = 1)
%     'precision'   sample precision (default = 'int16')
%    =========================================================================
%
% DEPENDENCIES:
%
%   LoadBinary
%   xxx


% Copyright (C) 2015 Adrien Peyrache, includes code from John D Long II
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.


function ConvertBasler2Whl(fbasename,digitalSig,varargin)


    pos = Process_DetectLED([fbasename '.avi']);
    dat = LoadBinary([fbasename '_digitalin.avi'],'nchannels',1,'channels',1);
    t = [0:length(dat)-1]/20000;
    
    df = diff(dat)<-1;
    frameT = t(df);
    
    if length(frameT)<size(pos,1);
        warning('Too many video frames!')
        pos = pos(1:length(frameT),:);
    elseif length(frameT)>size(pos,1);
        keyboard
    end
    
    
    
    
    
    

end

