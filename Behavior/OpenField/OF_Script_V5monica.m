%OF_Script MGP & LB March 2019
% This script will analyze open field experiment xy data. All measures are
% saved to a table and can be analyzed with AnalyzeOF.
%Input assumptions:
%
%
%Pipeline Description:
%
%   raw coordinates --> smoothed/transformed to cm --> calculate measures:
%
%       path length (cm)
%       instananeous velocity (cm/s)
%       dwell time for zones (s)
%       rate map (plotting & finding home base)
%       Home base coordinates
%       Number of stops
%       Circuitiy of path segments for stops
%       Proximity of stops to home base centroid
%
%   --> measures stored in table for later analysis.

%% Adding toolboxes, loading variables and setting frame rate (fr), ratemap bin size (binsize).

%MUST USE LABCOMPUTER 1 TO RUN THIS SCRIPT
com='D:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis\postprocess.m';
com=strsplit(com,filesep);
basedir=[com{1},filesep,'Users',filesep,com{3},filesep,'GoogleDrive',filesep,'MatlabDir'];
addpath([basedir,filesep,'chronux_2_11',filesep,'fly_track',filesep,'FAnalyze',filesep,'functions'],...
    [basedir,filesep,'BClarkToolbox',filesep,'Analysis'],...
    genpath([basedir,filesep,'BClarkToolbox',filesep,'Analysis',filesep,'Utils']),...
    [basedir,filesep,'BClarkToolbox',filesep,'Analysis',filesep,'Visualize']);

%SET THIS PATH TO WHERE YOUR PARAMS MAT FILE IS KEPT
load('d:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\PAE_OF\params_v1.mat') %%loads table containing the raw coords and xMax/Min and yMax/Min
clearvars -except params

fr=30; %frame rate
binsize=3; %to use for HB detection
upFactor=15;
params.HBcoords{1}=[];
params.dwellQuad{1}=[];
params.transcoords{1}=[];%initialize column in table for transformed xy coordinates
params.HBcenter{1}=[];
params.circuity{1}=[];
params.proximity{1}=[];

%% SMOOTH & TRANSFORM COORDINATES

%Fixing tracking errors with FixPos and transforming coordinates into cm,
% centered at 0,0

for i = 1:size(params.VideoName,1) %%size gives 2 values, # of rows and columns, by subtracting 1 we keep it on the interval 1:296
    
    params.ts{i}=linspace(0,size(params.rawcoords{i}(:,1),1)/fr, size(params.rawcoords{i}(:,1),1))'; %creates timestamps for Fixing tracking errors
    [trackData(:,1),trackData(:,2)] = FixPos(params.rawcoords{i}(:,1),params.rawcoords{i}(:,2),params.ts{i}); %smooths & reduces error from tracking
    [params.transcoords{i}(:,1),params.transcoords{i}(:,2)] = transformCoordinates(params.dia(i),params.Xmax(i),params.Xmin(i),params.Ymax(i),params.Ymin(i),trackData(:,1),trackData(:,2)); %Converts the coordinates into CM
    params.transcoords{i,1}(~any(~isnan(params.transcoords{i,1}), 2),:)=[];
    
    clear trackData
end

%% Compute Measures
param_idx=params.VideoName; %serves as index for getParam
for j=1:size(param_idx,1)
    
    %% Maze Parameters;
    
    % compute total distance and instantaneous velocity
    [ params.pathL{j},params.pathIV{j}] = compute_OFpathCalc(params.transcoords{j}(:,1),params.transcoords{j}(:,2),fr);
    
    %calculate verticies for quadrants
    quadrants=createZones([0,0],params.dia(j),'numQuad',16,'fig',0);
    
    %Calculate dwell time per quadrant
    for i=1:size(quadrants,1)-1
        tempXin=[0 quadrants(i,1) quadrants(i+1,1)];
        tempYin=[0 quadrants(i,2) quadrants(i+1,2)];
        [in,~]=inpolygon(params.transcoords{j,1}(:,1),params.transcoords{j,1}(:,2),tempXin,tempYin);
        params.dwellQuad{j}(1,i)=sum(in)/fr; %
        clear tempXin tempYin in
    end
    
    %Calculate verticies for wall area (%80)
    outsideWall=createZones([0,0],params.dia(j),'type','annulus','fig',0);
    
    %Calculate dwell time per annulus
    [in,~]=inpolygon(params.transcoords{j}(:,1),params.transcoords{j}(:,2),outsideWall(:,1),outsideWall(:,2));
    params.dwellOutsideWall{j}=sum(~in)/fr;
    
    %Find segments of path inbetween stops
    % FIND INDEX FOR VELOCITY FILTER (<3cm/second for 1 frame)
    params.stopIdx{j}=contiguousframes(params.pathIV{j}<3,30);
    [start,ends,~]=findgroups(params.stopIdx{j});
    
    %This finds coords for stopping
    for ii=1:length(start)
        motionless{ii}=[params.transcoords{j}(start(ii):ends(ii),1),params.transcoords{j}(start(ii):ends(ii),2)];
        %         plot(params.transcoords{j}(start(ii):ends(ii),1),params.transcoords{j}(start(ii):ends(ii),2),'*r');hold on
    end
    
    params.stops{j}=motionless;
  
    
    %This finds segments for running
    params.runIdx{j}=contiguousframes(params.pathIV{j}>3,15);
    [start,ends,~]=findgroups(params.runIdx{j});
    
    figure;
    for ii=1:length(start)
        motion{ii}=[params.transcoords{j}(start(ii):ends(ii),1),params.transcoords{j}(start(ii):ends(ii),2)];
        %         plot(params.transcoords{j}(start(ii):ends(ii),1),params.transcoords{j}(start(ii):ends(ii),2));hold on
    end
    
    
    %Compute circuity between segments
    for iii=1:length(motion)
        [params.circuity{j}(1,iii)]=circuity(motion{iii}(:,1),motion{iii}(:,2));
    end
    
    params.start{j}=motion;
    clear motion
    
    %Creates bin edges for heatmap
    xedge=linspace(-(params.dia(j)/2),(params.dia(j)/2),round(params.dia(j)/binsize));
    yedge=linspace(-(params.dia(j)/2),(params.dia(j)/2),round(params.dia(j)/binsize));
    
    %Bin coordinates for ratemap and apply guassian filter to smooth
    [map] = histcounts2(params.transcoords{j,1}(:,1),params.transcoords{j,1}(:,2),xedge,yedge);
    map=flipud(imrotate(map,90));
    map=map/fr;
    params.rateMap{j} = imgaussfilt(map, 1.5);  %STORE THIS!
    
    clear xedge yedge map
    
    %finds homebase coordinates
    [BW,maskedImage,x,y,fieldarea,params.upHBmap{j}] = segmentImage('map',params.rateMap{j},'upscalefac',upFactor); %Store area
    
    %rescales home base coordinates for application with transcoords
    for f=1:length(x)
        x_temp=rescale([x{f},upFactor+1,size(params.upHBmap{j},2)-upFactor],-(params.dia(j)/2),(params.dia(j)/2));
        x{f}=x_temp(1:end-2)';
        y_temp=rescale([y{f},upFactor+1,size(params.upHBmap{j},2)-upFactor],-(params.dia(j)/2),(params.dia(j)/2));
        y{f}=y_temp(1:end-2)';
    end
    
    clear start ends
    
    for f=1:size(x,1)
        params.HBcoords{j}{f,1}=x{f}; %save HB area coordinates for x
        params.HBcoords{j}{f,2}=y{f}; %save HB area coordinates for y
    end
    %Calculate centroid of HB
    for hb=1:size(params.HBcoords{j},1)
        [ ~,~, tempC] = min_encl_ellipsoid( params.HBcoords{j,1}{hb, 1},params.HBcoords{j,1}{hb, 2});
        params.HBcenter{j}(hb,1)=tempC(1,1);
        params.HBcenter{j}(hb,2)=tempC(2,1);
    end
    
    
    %%calculating proximity of HB center from the stop
    for num=1:length(motionless)
        x1=params.HBcenter{j}(1,1);
        x2=motionless{1,num}(1,1);
        y1=params.HBcenter{j}(1,2);
        y2=motionless{1,num}(1,2);
        
        params.proximity{j}(1,num)=sqrt((x1+x2)^2)-((y1+y2)^2);
    end
    
      clear motionless
      close all
      
    
    %
    %         %% save fig
    %         saveas(gcf,param_idx{j},'jpeg');
end

