%OF_preprocess MGP & LB March 2019; updated by LB May 2019
% This script will preprocess raw xy coordinates from open field
% experiments. 
%Input assumptions: 
%       A. Data are stored in a table named 'params' organized with the following columns:
%           1. PathName (string variable, path to video data)
%           2. VideoName (string variable - name of video file)
%           3. dia (diameter of open field used in experiment)
%           4 - 7. Xmax, Xmin, Ymax, Ymin (respectively). boundardy of open
%           field used to transform coordinates. 
%           8. rawcoords - nested cell array of Anymaze output.
%           9. structured cell array of xy coordinates for boundary of cue
%           used. 
%       B. params can be created by importing anymaze csv files. 
%
%Pipeline Description:
%
%   raw coordinates --> smoothed/transformed to cm --> calculate measures:
%


%% Adding toolboxes, loading variables and setting frame rate (fr), ratemap bin size (binsize).
%Navigate to folder containing ephys_tools
s=what('ephys_tools');
cd(s.path)

%Initialize Paths 
addpath(fullfile('Behavior','OpenField','MeasureFunctions'),...
    fullfile('Behavior','OpenField'),...
    fullfile('Behavior'),...
    fullfile('Utils'),...
    fullfile('Visualize'),...
    fullfile('external_packages','chronux_2_11'),...
    fullfile('external_packages','Colormaps'))

cd('D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data')
load('params_V11.mat') %%loads table containing the raw coords and xMax/Min and yMax/Min
clearvars -except params

fr=30; %frame rate
binsize=3; %to use for HB detection
upFactor=15; %used in segmentimage - this upscales the heatmap to get a smooth boundary for home base. 

%Initialze table for measures 
params.HBcoords{1}=[];
params.dwellQuad{1}=[];
params.transcoords{1}=[];%initialize column in table for transformed xy coordinates
params.HBcenter{1}=[];
params.circuity{1}=[];
params.proximity{1}=[];
params.transcue{1}=[];

%% SMOOTH & TRANSFORM COORDINATES

%Fixing tracking errors with FixPos and transforming coordinates into cm,
% centered at 0,0
for i = 1:size(params.PathName,1)-1 %%size gives 2 values, # of rows and columns, by subtracting 1 we keep it on the interval 1:296
    params.ts{i}=linspace(0,size(params.rawcoords{i}(:,1),1)/fr, size(params.rawcoords{i}(:,1),1))'; %creates timestamps for Fixing tracking errors
    [trackData(:,1),trackData(:,2)] = FixPos(params.rawcoords{i}(:,1),params.rawcoords{i}(:,2),params.ts{i}); %smooths & reduces error from tracking
    [params.transcoords{i}(:,1),params.transcoords{i}(:,2)] = transformCoordinates(params.dia(i),params.Xmax(i),params.Xmin(i),params.Ymax(i),params.Ymin(i),trackData(:,1),trackData(:,2)); %Converts the coordinates into CM
    params.transcoords{i,1}(~any(~isnan(params.transcoords{i,1}), 2),:)=[];
    clear trackData
end

%%Transform Cue Coords
for j=1:size(params.PathName,1)

    fields = fieldnames(params.cue{j});
    it=1;
    for i = 1:2:(size(fields,1)) %%size gives 2 values, # of rows and columns, by subtracting 1 we keep it on the interval 1:296
        [tempCueX(it,1),tempCueY(it,1)] =...
            transformCoordinates(params.dia(j),params.Xmax(j),params.Xmin(j),...
            params.Ymax(j),params.Ymin(j),params.cue{j}.(fields{i}){:},params.cue{j}.(fields{i+1}){:}); %Converts the coordinates into CM
        it=it+1;
    end
    params.transcue{j}(:,1)=tempCueX;
    params.transcue{j}(:,2)=tempCueY;
    
    clear tempCueX tempCueY
end

