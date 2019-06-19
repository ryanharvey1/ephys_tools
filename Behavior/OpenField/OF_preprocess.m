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
addpath(fullfile('external_packages','chronux_2_11'),...
    fullfile('external_packages','Colormaps'))

cd('D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data')
load('params_V14.mat') %%loads table containing the raw coords and xMax/Min and yMax/Min
clearvars -except params

fr=30; %frame rate

%Initialze table for measures
params.HBcoords{1}=[];
params.dwellQuad{1}=[];
params.noseCM{1}=[];
params.headCM{1}=[];
params.backCM{1}=[];
params.buttCM{1}=[];%initialize column in table for transformed xy coordinates
params.HBcenter{1}=[];
params.circuity{1}=[];
params.proximity{1}=[];
params.cueCM{1}=[];
params.hbOcc{1}=[];
params.hbVel{1}=[];
params.hbOcc{1}=[];
params.hbVel{1}=[];
params.hbXYidx{1}=[];
params.entries{1}=[];
params.HBclass{1}=[];
params.HBBound{1}=[];
params.slowInHB{1}=[];

%% SMOOTH & TRANSFORM COORDINATES

%transforming coordinates into cm, centered at 0,0
for i = 1:size(params.subID,1)
    
    [params.noseCM{i}(:,1),params.noseCM{i}(:,2)] = transformCoordinates...
        (params.dia{i},params.Xmax{i},params.Xmin{i},params.Ymax{i},params.Ymin{i},...
        params.nose{i}(:,1),params.nose{i}(:,2)); %Converts the coordinates into CM
    
    [params.headCM{i}(:,1),params.headCM{i}(:,2)] = transformCoordinates...
        (params.dia{i},params.Xmax{i},params.Xmin{i},params.Ymax{i},params.Ymin{i},...
        params.head{i}(:,1),params.head{i}(:,2)); %Converts the coordinates into CM
    
    [params.backCM{i}(:,1),params.backCM{i}(:,2)] = transformCoordinates...
        (params.dia{i},params.Xmax{i},params.Xmin{i},params.Ymax{i},params.Ymin{i},...
        params.back{i}(:,1),params.back{i}(:,2)); %Converts the coordinates into CM
    
    [params.buttCM{i}(:,1),params.buttCM{i}(:,2)] = transformCoordinates...
        (params.dia{i},params.Xmax{i},params.Xmin{i},params.Ymax{i},params.Ymin{i},...
        params.butt{i}(:,1),params.butt{i}(:,2)); %Converts the coordinates into CM
    
    clear trackData
end

%%Transform Cue Coords
for j=1:size(params.subID,1)
    
    if isempty(params.cueCoords{j})
        params.cueCM{j}(:,1)=NaN;
        params.cueCM{j}(:,2)=NaN;
        continue
    end
  
    [params.cueCM{j}(:,1),params.cueCM{j}(:,2)] =...
        transformCoordinates(params.dia{j},params.Xmax{j},params.Xmin{j},...
        params.Ymax{j},params.Ymin{j},params.cueCoords{j}(:,1),params.cueCoords{j}(:,2)); %Converts the coordinates into CM
end



