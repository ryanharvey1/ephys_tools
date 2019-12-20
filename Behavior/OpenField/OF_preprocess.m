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
%   raw coordinates --> transformed to cm --> calculate measures:
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
params.HBcoords{1}=[]; %home base coordinates
params.dwellQuad{1}=[]; %time spent in maze qudrants 
params.noseCM{1}=[]; %xy coordinates for nose in centimeters
params.headCM{1}=[]; %xy coordinates for head in centimeters
params.backCM{1}=[]; %xy coordinates for back in centimeters
params.buttCM{1}=[]; %xy coordinates for butt in centimeters
params.HBcenter{1}=[]; %xy coordiantes for center of home base(s)
params.circuity{1}=[]; %circuity of path segments when animal is moving
params.proximity{1}=[]; %proximity of 
params.cueCM{1}=[]; %coordinates for the cue in centimeters
params.hbOcc{1}=[]; %time spent in home base
params.hbVel{1}=[]; %average velocity when rat is in home base
params.hbXYidx{1}=[]; %index for when animal is in a home base
params.entries{1}=[]; % number of times animal enters the home base (must be there for a min of 4 seconds)
params.HBclass{1}=[]; % ratio between time spent in hb and mean velocity in hb. For hb classification. 
params.HBBound{1}=[]; % xy coordintes of the hb boundary
params.slowInHB{1}=[]; %time spent slow in home base
params.HBstops{1}=[]; % number of stops in the home base boundary (determined by each rat's home base)
params.HBdist2Cue{1}=[]; %distance from hb center to the cue
params.distHBstop{1}=[]; %average distance from hb center to stopping location
params.closeHBstop{1}=[]; %Number of stops within a 25cm diameter from HB center
params.pathLQuad{1}=[]; %path length of distance in pie shaped quadrant
params.pathDist{1}=[]; %distance between consecutive points in path. 
params.pathIVQuad{1}=[]; %instantaneous linear velocity within quadrant. 
params.angVelQuad{1}=[]; %instantaneous angular velocity within quadrant
params.stopQuad{1}=[]; %time stopped in quadrant. 
params.numstopQuad{1}=[]; %number of times animal stopped in quadrant.
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



