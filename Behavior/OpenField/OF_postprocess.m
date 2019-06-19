%OF_postprocessing

%Initialize Paths 
addpath(fullfile('Behavior','OpenField','MeasureFunctions'),...
    fullfile('Behavior','OpenField'),...
    fullfile('Behavior'),...
    fullfile('Utils'),...
    fullfile('Visualize'),...
    fullfile('external_packages','chronux_2_11'),...
    fullfile('external_packages','Colormaps'))


cd('D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data')
load('params_V15.mat') %%loads table created from OF_preprocess. 
clearvars -except params


%% _________________________Compute Measures_______________________________

% Measures using 'Back' coordinates: 
% pathL: path length (cm)
% pathIV: instantaneous veloctiy (cm/s)
% dwellQuad: time in equally spaced 16 pie shaped segments (s)
% dwellOutsideWall: time spent in outer 20% of maze (s)
% stopIdx: saves indicies of when the animals movement speed drops below
% 3cm/s for at least one second.
% 
% Measures using 'Nose' & 'Head' coordinates: 
%
% HD: allocentric head direction 
% cueInspection: Time spent with the rats nose or head within <5cm from the
% cue location. 
% 
fr=30; %frame rate
binsize=3; %to use for HB detection
upFactor=15; %used in segmentimage - this upscales the heatmap to get a smooth boundary f

param_idx=params.subID; %serves as index for getParam
params.hbOcc{1}=[];
params.hbVel{1}=[];
params.hbXYidx{1}=[];
params.entries{1}=[];
params.HBclass{1}=[];
params.HBBound{1}=[];
params.slowInHB{1}=[];
%Loop through subjects
for j=1:size(param_idx,1)
    
    %Grab variables for easier use
    back=params.backCM{j};
    nose=params.noseCM{j};
    head=params.headCM{j};
    butt=params.buttCM{j};

    % compute total distance (pathL) and instantaneous velocity (pathIV)
    [ params.pathL{j},params.pathIV{j}] = compute_OFpathCalc(back(:,1),back(:,2),fr);
    
    params.runSpeed{j}=nanmean(params.pathIV{j}(params.pathIV{j}>=3));
    
    
    %Head direction 
    params.HD{j}=wrapTo360(fixNLXangle(rad2deg(atan2(head(:,2)-nose(:,2),...
        head(:,1)-nose(:,1))),round(0.1667*30)));

    %calculate verticies for quadrants
    quadrants=createZones([0,0],params.dia{j},'numQuad',16,'fig',0); %16 pie shaped segments
    
    %Calculate dwell time per quadrant
    for i=1:size(quadrants,1)-1
        tempXin=[0 quadrants(i,1) quadrants(i+1,1)];
        tempYin=[0 quadrants(i,2) quadrants(i+1,2)];
        [in,~]=inpolygon(back(:,1),back(:,2),tempXin,tempYin);
        params.dwellQuad{j}(1,i)=sum(in)/fr; %
        clear tempXin tempYin in
    end
    
    %Calculate verticies for wall area (%80)
    outsideWall=createZones([0,0],params.dia{j},'type','annulus','fig',0); %default annulus size is 80% for createZones
    
    %Calculate dwell time per annulus
    [in,~]=inpolygon(back(:,1),back(:,2),outsideWall(:,1),outsideWall(:,2));
    params.dwellOutsideWall{j}=sum(~in)/fr;
    
    %Find segments of path inbetween stops
    % FIND INDEX FOR VELOCITY FILTER (<3cm/second for 1 second)
    params.stopIdx{j}=contiguousframes(params.pathIV{j}<3,30);
    [start,ends,~]=findgroups(params.stopIdx{j});
    
    %This finds coords for stopping
    for ii=1:length(start)
        motionless{ii}=[back(start(ii):ends(ii),1),back(start(ii):ends(ii),2)];
    end
    
    params.stops{j}=motionless; clear motionless 
    
    %%Create Number of Stops
    params.NumStops{j}= size(params.stops{j},2);
    
    %This finds segments for running
    params.runIdx{j}=contiguousframes(params.pathIV{j}>3,15);
    [start,ends,~]=findgroups(params.runIdx{j});
    
    figure;
    for ii=1:length(start)
        motion{ii}=[back(start(ii):ends(ii),1),back(start(ii):ends(ii),2)];
        [ segPL{ii},segPathIV{ii}] = compute_OFpathCalc(back(start(ii):ends(ii),1),...
            back(start(ii):ends(ii),2),fr);
        %         plot(params.transcoords{j}(start(ii):ends(ii),1),params.transcoords{j}(start(ii):ends(ii),2));hold on
    end
    params.runs{j}=motion;   clear motion
    params.segPL{j}=segPL; clear segPL
    params.segPathIV{j}=segPathIV; clear segPathIV
     
    clear ends start
    
    %Compute circuity between segments
    for iii=1:size(params.runs{j},2)
        [params.circuity{j}(1,iii)]=circuity(params.runs{j}{iii}(:,1),params.runs{j}{iii}(:,2));
    end
    
    params.meanCirc{j}=nanmean(params.circuity{j}(1,:)); %compute average circuity
    
    
    params.NumRuns{j}= size(params.runs{j},2);
    
    %Creates bin edges for heatmap
    xedge=linspace(-(params.dia{j}/2),(params.dia{j}/2),round(params.dia{j}/binsize));
    yedge=linspace(-(params.dia{j}/2),(params.dia{j}/2),round(params.dia{j}/binsize));
    
    %Bin coordinates for heat map and apply guassian filter to smooth
    [map] = histcounts2(back(:,1),back(:,2),xedge,yedge);
    map=flipud(imrotate(map,90));
    map=map/fr;
    params.rateMap{j} = imgaussfilt(map, 1.5);  %STORE THIS!
    
    clear xedge yedge
    
    %Use meshgrid to serve as basis for logical mask. 
    imageSizeX = size(map,1);
    imageSizeY = size(map,2);
    [columnsInImage,rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
    
    clear imageSizeX imageSizeY
    
    % Next create the circle in the image.
    centerX = median(1:size(map,1)); centerY = median(1:size(map,2)); radius = median(1:size(map,2));
    circlePixels = (rowsInImage - centerY).^2 ...
        + (columnsInImage - centerX).^2 <= radius.^2;
    
    clear rowsInImage columnsInImage 
    
    map(~circlePixels)=NaN; %indicate area outside the maze by labelling with NaN
    params.searchArea{j}=sum(sum(map>0))/sum(sum(~isnan(map))); %Calculate the proportion of bins occupied by animal.
    
    clear xedge yedge map circlePixels
    
    %finds high occupancy coordinates
    [~,~,x,y,params.fieldarea{j},params.upHBmap{j}] = segmentImage('map',params.rateMap{j},'upscalefac',upFactor); %Store area
    
  
    %get home base measures
    for f=1:length(x)
        %rescale coordinates back to pool size
        x_temp=rescale([x{f},upFactor+1,size(params.upHBmap{j},2)-upFactor],-(params.dia{j}/2),(params.dia{j}/2));
        y_temp=rescale([y{f},upFactor+1,size(params.upHBmap{j},2)-upFactor],-(params.dia{j}/2),(params.dia{j}/2));
        tempIn=inpolygon(params.backCM{j,1}(:,1),params.backCM{j,1}(:,2),x_temp(1:end-2)',y_temp(1:end-2)');    
        tempOut=contiguousframes(tempIn,60); %has to be inside of hb for at least 2 sec to count as entry
        [~,~,params.entries{j}{1,f}]=findgroups(tempOut);
        params.hbOcc{j}(1,f)=nansum(tempIn)/fr; %total time in home base
        params.hbVel{j}(1,f)=nanmean(params.pathIV{j,1}(tempIn(1:end-1,1),1)); %remove last tempIn idx to accomodate PathIV length
        params.hbXYidx{j}{1,f}=tempIn;
        params.HBBound{j}{1,f}=[x_temp(1:end-2)',y_temp(1:end-2)'];
        params.HBcoords{j}{1,f}=[params.backCM{j,1}(tempIn,1),params.backCM{j,1}(tempIn,2)];
        params.slowInHB{j}{1,f}=nansum(tempIn(1:end-1)& params.stopIdx{j})/fr; %time being slow in homebase
        params.HBclass{j}{1,f}=params.slowInHB{j}{1,f}/params.hbOcc{j}(1,f); %proportion of time being slow in hb
    end
  
     clear x_temp y_temp tempIn tempOut x y
     
    %Calculate centroid of high occupancy coordinates
    for hb=1:size(params.HBBound{j},2)
        [ ~,~, tempC] = min_encl_ellipsoid( params.HBBound{j}{1, hb}(:,1),params.HBBound{j}{1, hb}(:,2));
        params.HBcenter{j}{1,hb}=[tempC(1,1),tempC(2,1)]; %x y coords
    end
%     
    %Calculate distance measures between high occupancy coordinates centers
    if size(params.HBcenter{j},1)>1
        params.avgHBdist{j}=mean(pdist(params.HBcenter{j}));
        params.maxHBdist{j}=max(pdist(params.HBcenter{j}));
        params.minHBdist{j}=min(pdist(params.HBcenter{j}));
    else
        params.avgHBdist{j}=NaN;
        params.maxHBdist{j}=NaN;
        params.minHBdist{j}=NaN;
    end
    
    %%calculating proximity of high occupancy coordinates center from the stop
%  
%     
%     %Find number of stops with 25cm diameter of hb center
%     params.stopsNearHB{j}=sum(params.proximity{j}<25);
% 
%     %Find proximity of home base to cue edge
%     k = convhull(params.cueCM{j}(:,1),params.cueCM{j}(:,2));
%     cueBoundary=[params.cueCM{j}(k,1),params.cueCM{j}(k,2)];
%     
%     distances = sqrt(sum(bsxfun(@minus, cueBoundary, [xinzone(j),yinzone(j)]).^2,2));
%     closestwall = boundary(distances==min(distances),:);
%     
close all
end
