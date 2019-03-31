function [ track_length ] = TrackLength( path )
% TrackLength Obtains tracklength based on certain cut off dates
% Ryan Harvey 3/22/2017
%
% SPECIFY TRACK LENGTH >>> 90 (pre 7/14/16) or 120 (7/14/16 and on)
year=strsplit(path,filesep);
year=strsplit(char(year(end)),'-');
if strcmp(year(1),num2str(2016)) || strcmp(year(1),num2str(2015))
    SplitPath=strsplit(path,'-'); % SPLIT PATH UP BY -
    SplitPath2=strsplit(char(SplitPath(3)),'_'); % SPLIT PATH UP BY _
    if str2double(SplitPath(2))<=7; % IF MONTH IS 7 OR LESS
        track_length = 90;
        if str2double(SplitPath(2))==7;
            if str2double(SplitPath2(1))<13; % IF DAY IS LESS THAN 14
                track_length = 90;
            else str2double(SplitPath2(1))>13;
                track_length = 120;
            end
        end
    else str2double(SplitPath(2))>7;
        track_length = 120;
    end
else
    track_length = 120;
end
end

