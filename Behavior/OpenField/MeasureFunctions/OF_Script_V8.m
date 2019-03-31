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
com='F:\Users\BClarkLab\GoogleDrive\MatlabDir';
com=strsplit(com,filesep);
basedir=[com{1},filesep,'Users',filesep,com{3},filesep,'GoogleDrive',filesep,'MatlabDir'];
addpath([basedir,filesep,'chronux_2_11',filesep,'fly_track',filesep,'FAnalyze',filesep,'functions'],...
    [basedir,filesep,'BClarkToolbox',filesep,'Analysis'],...
    genpath([basedir,filesep,'BClarkToolbox',filesep,'Analysis',filesep,'Utils']),...
    [basedir,filesep,'BClarkToolbox',filesep,'Analysis',filesep,'Visualize']);

load('params_V6.mat') %%loads table containing the raw coords and xMax/Min and yMax/Min
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
    if  ~contains(params.PathName{j},"day1_lgOF")
        
        continue
    end
    
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

%% Compute Measures
param_idx=params.PathName; %serves as index for getParam
for j=1:size(param_idx,1)
    
    if ~contains(param_idx{j},'lgOF') % to bypass small open field for now LB March 2019
        continue
    end
    
    close all
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
    
    params.stops{j}=motionless; clear motionless 
    
    %%Create Number of Stops
    params.NumStops{j}= size(params.stops{j},2);
    
    %This finds segments for running
    params.runIdx{j}=contiguousframes(params.pathIV{j}>3,15);
    [start,ends,~]=findgroups(params.runIdx{j});
    
    figure;
    for ii=1:length(start)
        motion{ii}=[params.transcoords{j}(start(ii):ends(ii),1),params.transcoords{j}(start(ii):ends(ii),2)];
        %         plot(params.transcoords{j}(start(ii):ends(ii),1),params.transcoords{j}(start(ii):ends(ii),2));hold on
    end
    params.runs{j}=motion;   clear motion
    
    %Compute circuity between segments
    for iii=1:size(params.runs{j},2)
        [params.circuity{j}(1,iii)]=circuity(params.runs{j}{iii}(:,1),params.runs{j}{iii}(:,2));
    end
    
    params.meanCirc{j}=nanmean(params.circuity{j}(1,:)); %compute average circuity
    
    
    params.NumRuns{j}= size(params.runs{j},2);
    
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
    
    params.HBcoords{j}(:,1)=x; %save HB area coordinates for x
    params.HBcoords{j}(:,2)=y; %save HB area coordinates for y
    
    %Calculate centroid of HB
    for hb=1:size(params.HBcoords{j},1)
        [ ~,~, tempC] = min_encl_ellipsoid( params.HBcoords{j,1}{hb, 1},params.HBcoords{j,1}{hb, 2});
        params.HBcenter{j}(hb,1)=tempC(1,1);
        params.HBcenter{j}(hb,2)=tempC(2,1);
    end
    
    
    %%calculating proximity of HB center from the stop
    for q=1:size(params.HBcoords{j},1)
        max_HB_idx=length(params.HBcoords{j}{q});
    end
    [~,max_HB_idx]=max(max_HB_idx);
    
    for num=1:size(params.stops{j},2)
        params.proximity{j}(1,num)=sqrt((params.HBcenter{j}(max_HB_idx,1)-params.stops{j}{1,num}(1,1))^2+(params.HBcenter{j}(max_HB_idx,2)-params.stops{j}{1,num}(1,2))^2);
    end
    
    %Find number of stops with 20cm diameter of hb center
    params.stopsNearHB{j}=sum(params.proximity{j}<20);
    
    
end


%Code Graveyard

%RESCALES XY COORDINATES TO MATCH THE UPSCALED RATEMAP OBTAINED IN
%SEGMENTIMAGE.
%         rescale_x=rescale([params.transcoords{j,1}(:,1);-(params.dia(j)/2);(params.dia(j)/2)],upFactor+1,size(X,2)-upFactor);
%         rescale_y=rescale([params.transcoords{j,1}(:,2);-(params.dia(j)/2);(params.dia(j)/2)],upFactor+1,size(X,2)-upFactor);
%
%         rescale_x=rescale_x(1:end-2);
%         rescale_y=rescale_y(1:end-2);
%         figure;imagesc(X);hold on;plot(rescale_x,rescale_y,'k')
%
%         rescale_y=rescale([params.transcoords{j,1}(:,2);-(params.dia(j)/2);(params.dia(j)/2)],upFactor+1,size(params.upHBmap{j},2)-upFactor);
%
%         rescale_x=rescale_x(1:end-2);
%         rescale_y=rescale_y(1:end-2);
%         figure;imagesc(X);hold on;plot(rescale_x,rescale_y,'k')
%
%
%
%
