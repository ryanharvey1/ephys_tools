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
load('params_V17.mat') %%loads table created from OF_preprocess.
clearvars -except params
params.time2HB{1}=[];

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
%Loop through subjects
for j=1:size(param_idx,1)
    
    %Grab variables for easier use
    back=params.backCM{j};
    nose=params.noseCM{j};
    head=params.headCM{j};
    butt=params.buttCM{j};
    ts=params.ts{j};
    
    % compute total distance (pathL) and instantaneous velocity (pathIV)
    [ params.pathL{j},params.pathIV{j},params.pathDist{j}] = compute_OFpathCalc(back(:,1),back(:,2),fr);
    params.runSpeed{j}=nanmean(params.pathIV{j}(params.pathIV{j}>=3));
%     
%     Find segments of path inbetween stops

%     FIND INDEX FOR VELOCITY FILTER (<3cm/second for 1 second)
    params.stopIdx{j}=contiguousframes(params.pathIV{j}<3,30);
    [startStop,endStop,~]=findgroups(params.stopIdx{j});
    
%     This finds coords for stopping
    for ii=1:length(startStop)
        motionless{ii}=[back(startStop(ii):endStop(ii),1),back(startStop(ii):endStop(ii),2)];
        timeMotionless{ii}=size(motionless{ii},1)/fr;
        tsStop{ii}=ts(startStop(ii):endStop(ii),1);
    end
    
    params.stops{j}=motionless; clear motionless
    params.timeStopped{j}=timeMotionless; clear timeMotionless
    params.tsStop{j}=tsStop; clear tsStop
    
%     Create Number of Stops
    params.NumStops{j}= size(params.stops{j},2);
    
    %Head direction
    params.HD{j}=wrapTo360(fixNLXangle(rad2deg(atan2(head(:,2)-nose(:,2),...
        head(:,1)-nose(:,1))),round(0.1667*30)));
    
    params.angVel{j}=insta_angvel(params.HD{j}',fr);
    
    %calculate verticies for quadrants
    quadrants=createZones([0,0],params.dia{j},'numQuad',16,'fig',0); %16 pie shaped segments
    
    %Find center of mass for stops.
    for ii=1:size(params.stops{j},2)
        temp=[params.stops{j}{1,ii}(:,1),params.stops{j}{1,ii}(:,2)];
        temp(any(isnan(temp),2),:)=[];
        if ~isempty(temp) && size(temp,1)>1
            [ ~,~, stopCenter] = min_encl_ellipsoid(temp(:,1),temp(:,2));
        else
            stopCenter(1,1)=NaN;
            stopCenter(2,1)=NaN;
        end
        stop(ii,1)=stopCenter(1,1);
        stop(ii,2)=stopCenter(2,1);
    end
    
    params.stopCenter{j}=stop; 
    clear stop
    
    %Calculate dwell time per quadrant
    for i=1:size(quadrants,1)-1
        tempXin=[0 quadrants(i,1) quadrants(i+1,1)];
        tempYin=[0 quadrants(i,2) quadrants(i+1,2)];
        [in,~]=inpolygon(back(:,1),back(:,2),tempXin,tempYin);
        params.dwellQuad{j}(1,i)=sum(in)/fr; %
        params.pathLQuad{j}(1,i)= sum(params.pathDist{j}(in(1:end-1)& params.pathIV{j}>3 ,1)); %mean distance when animal is in quadrant and moving >3cm/s
        params.pathIVQuad{j}(1,i)= nanmean(params.pathIV{j}(in(1:end-1),1)); %mean linear velocity
        params.numstopQuad{j}(1,i)=sum(inpolygon(params.stopCenter{j}(:,1),params.stopCenter{j}(:,2),tempXin,tempYin));
        params.stopQuad{j}(1,i)=sum(in(1:end-1,:) & params.stopIdx{j})/fr; %stop duration
        params.angVelQuad{j}(1,i)=nanmean(abs(params.angVel{j}(in(1:end-1,:))));
        clear tempXin tempYin in
    end
    
    %Create annuli for center and outter edge of maze
    outsideDwell=createZones([0,0],params.dia{j},'type','annulus','fig',0,'annulusSize',.6); %default annulus size is 80% for createZones
    centerDwell=createZones([0,0],params.dia{j},'type','annulus','fig',0,'annulusSize',.6); %default annulus size is 80% for createZones
    
    %Calculate dwell time for outter edge
    [in,~]=inpolygon(back(:,1),back(:,2),outsideDwell(:,1),outsideDwell(:,2));
    params.dwellOutsideWall{j}=sum(~in)/fr;
    
    clear in outsideDwell
    
    %Calculate dwell time for center maze
    [in,~]=inpolygon(back(:,1),back(:,2),centerDwell(:,1),centerDwell(:,2));
    params.centerDwell{j}=sum(in)/fr;
    
    clear in centerDwell
    
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
        params.slowInHB{j}{1,f}=nansum(tempIn(1:end-1)& params.stopIdx{j})/fr; % time being slow in homebase
        params.HBclass{j}{1,f}=params.slowInHB{j}{1,f}/params.hbOcc{j}(1,f); % proportion of time being slow in hb
        params.HBstops{j}{1,f}=nansum(inpolygon(params.backCM{j,1}(startStop,1),...
        params.backCM{j,1}(startStop,2),x_temp(1:end-2)',y_temp(1:end-2)')); % Find number of times animals initiated a start in the home base
        time2HB_idx=inpolygon(params.backCM{j,1}(startStop,1),...
        params.backCM{j,1}(startStop,2),x_temp(1:end-2)',y_temp(1:end-2)');
        firstIdx=find(time2HB_idx);
        time2HB_time=params.tsStop{j,1}{1,firstIdx(1)};
        params.time2HB{j}{1,f}=time2HB_time(1,1);
    end
    
    
    clear x_temp y_temp tempIn tempOut x y distances
    
    %Calculate centroid of high occupancy coordinates
    for hb=1:size(params.HBBound{j},2)
        [ ~,~, tempC] = min_encl_ellipsoid( params.HBBound{j}{1, hb}(:,1),params.HBBound{j}{1, hb}(:,2));
        params.HBcenter{j}{1,hb}=[tempC(1,1),tempC(2,1)]; %x y coords
        
    end
    
    %Average Proximity of stops from hb center
    for hb=1:size(params.HBcenter{j},2)
        tempHB=params.HBcenter{j}{1,hb};
        for i=1:size(params.stops{j},2)
            temp=params.stops{j}{1,i};
            dist=sqrt((tempHB(1,1)-temp(:,1)).^2+(tempHB(1,2)-temp(:,2)).^2);
            disttemp(i,1)=nanmean(dist);
            
        end
        params.distHBstop{j}{1,hb}=nanmean(disttemp); %average distance from stops 
        params.closeHBstop{j}{1,hb}=sum(disttemp<25); %number of stops within 25cm from hb center
        params.distHBstopAll{j}=disttemp; clear disttemp tempHB dist temp
    end
    
    %Calculate distance measures between high occupancy coordinates centers
    temp=[];
    if size(params.HBcenter{j},2)>1
        for hb=1:size(params.HBcenter{j},2)
            temp=[temp; params.HBcenter{j}{1,hb}];
        end
        params.avgHBdist{j}=nanmean(pdist(temp));
        params.maxHBdist{j}=max(pdist(temp));
        params.minHBdist{j}=min(pdist(temp));
    else
        params.avgHBdist{j}=NaN;
        params.maxHBdist{j}=NaN;
        params.minHBdist{j}=NaN;
    end
    clear temp
    %%calculating proximity of high occupancy coordinates center from the
    %%cue boundary
    if ~isempty(params.cueCoords{j})
        k = convhull(params.cueCM{j}(:,1),params.cueCM{j}(:,2));
        cueBoundary=[params.cueCM{j}(k,1),params.cueCM{j}(k,2)];
        
        for r=1:size(params.HBcenter{j},2)
            distances = sqrt(sum(bsxfun(@minus, cueBoundary, [params.HBcenter{j}{1,r}(:,1),params.HBcenter{j}{1,r}(:,2)]).^2,2));
            params.HBdist2Cue{j}{1,r} = unique(distances(distances==min(distances))); %Find minimum distance from cue boundary to hb center
        end
    else
         k = convhull(params.cueCM{j-1}(:,1),params.cueCM{j-1}(:,2));
        cueBoundary=[params.cueCM{j-1}(k,1),params.cueCM{j-1}(k,2)];
        
        for r=1:size(params.HBcenter{j},2)
            distances = sqrt(sum(bsxfun(@minus, cueBoundary, [params.HBcenter{j}{1,r}(:,1),params.HBcenter{j}{1,r}(:,2)]).^2,2));
            params.HBdist2Cue{j}{1,r} = unique(distances(distances==min(distances))); %Find minimum distance from cue boundary to hb center
        end
    end
    
    % Excursion Analysis 
   
    %This finds stops that occur in the home base boundary
    excursion{1}=[];
    params.stopIdx{j}=contiguousframes(params.pathIV{j}<3,30);
    [startStop,endStop,~]=findgroups(params.stopIdx{j});
    
    for hb=1:size(params.HBBound{j},2)
        endIn{hb} = inpolygon(back(endStop',1),back(endStop',2), params.HBBound{j}{1, hb}(:,1),params.HBBound{j}{1, hb}(:,2));
        startIn{hb} = inpolygon(back(startStop',1),back(startStop',2), params.HBBound{j}{1, hb}(:,1),params.HBBound{j}{1, hb}(:,2));
        
        for runs=1:size(params.runs{j},2)
            runIn = inpolygon(params.runs{j}{1,runs}(:,1),params.runs{j}{1,runs}(:,2), params.HBBound{j}{1, hb}(:,1),params.HBBound{j}{1, hb}(:,2));
            [row,~]=find(runIn,1,'first');
            if isempty(row)
                excursion{hb}(runs,1)=NaN;
                excursion{hb}(runs,2)=NaN;
                continue
            end
            excursion{hb}(runs,1)=params.runs{j}{1,runs}(row,1);
            excursion{hb}(runs,2)=params.runs{j}{1,runs}(row,2);
       
        end
     
    end
    
    params.startInHB{j}=endIn; clear endIn
    params.endInHB{j}=startIn; clear startIn
    %find time when moving greater than or equal to 3cm/s
    params.stopIdx{j}=contiguousframes(params.pathIV{j}>=3,30);
    [startMotion,endMotion,~]=findgroups(params.stopIdx{j});
    
    for ii=1:length(startMotion)
        motion{ii}=[back(startMotion(ii):endMotion(ii),1),back(startMotion(ii):endMotion(ii),2)];
        timeMotion{ii}=size(motion{ii},1)/fr;
        [ segPL,segPathIV,~] = compute_OFpathCalc(back(:,1),back(:,2),fr);
        tsMotion{ii}=ts(startMotion(ii):endMotion(ii),1);
    end
    
    
    params.runs{j}=motion;   clear motion
    params.timeMoving{j}=timeMotion; clear timeMotion
    params.segPL{j}=segPL; clear segPL
    params.segPathIV{j}=segPathIV; clear segPathIV
    params.tsMotion{j}=tsMotion; clear tsMotion
        
    %Compute circuity between segments
    for iii=1:size(params.runs{j},2)
        [params.circuity{j}(1,iii)]=circuity(params.runs{j}{iii}(:,1),params.runs{j}{iii}(:,2));
    end
    
    params.meanCirc{j}=nanmean(params.circuity{j}(1,:)); %compute average circuity
    params.NumRuns{j}= size(params.runs{j},2);
    
end
