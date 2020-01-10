%OF_postprocessing

%Initialize Paths
addpath(fullfile('Behavior','OpenField','MeasureFunctions'),...
    fullfile('Behavior','OpenField'),...
    fullfile('Behavior'),...
    fullfile('Utils'),...
    fullfile('Visualize'),...
    fullfile('external_packages','chronux_2_11'),...
    fullfile('external_packages','Colormaps'))


% cd('D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data')
% load('params_V18.mat') %%loads table created from OF_preprocess.
clearvars -except params

%% _________________________Compute Measures_______________________________

% Measures using 'Back' coordinates:
% pathL: path length (cm)
% pathIV: instantaneous veloctiy (cm/s)
% dwellQuad: time in equally spaced 16 pie shaped segments (s)
% dwellOutsideWall: time spent in outer 20% of maze (s)
% stopIdx: saves indicies of when the animals movement speed drops below
% 3cm/s for at least one second.
% runIdx: save indicies of running movements
%
% Measures using 'Nose' & 'Head' coordinates:
%
% HD: allocentric head direction
% cueInspection: Time spent with the rats nose or head within < 5 cm from the
% cue location.
%

%% Enter Setting 
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
    ts=params.ts{j};
 %% Basic locomotor Measures (path length, velocity, number of stops/stop duration, time spent in outer/inner part of maze)
 
    % compute total distance (pathL) and instantaneous velocity (pathIV)
    [params.pathL{j},params.pathIV{j},params.runSpeed{j},params.pathAC{j},params.mean_accel{j},params.pathDist{j}] = OF.path_measures(j,params,fr);
    
    %Head direction
    params.HD{j}=wrapTo360(fixNLXangle(rad2deg(atan2(head(:,2)-nose(:,2),...
        head(:,1)-nose(:,1))),round(0.1667*30)));
    params.angVel{j}=insta_angvel(params.HD{j}',fr);
    
    % Compute Stopping behaviors
    stop_measures = OF.stops(j,params,fr);
    
    params.stopIdx{j} = stop_measures.stopIdx; 
    params.stops{j} = stop_measures.stops;
    params.timeStopped{j} = stop_measures.timeStopped;
    params.tsStop{j} = stop_measures.tsStop;
    params.NumStops{j} = stop_measures.NumStops;
    params.stopCenter{j} = stop_measures.stopCenter; 
    
    % Quadrant Analysis - quantify OF measures given a quadrant 
    
    quad_measures = OF.quadrant(j,params,fr,16); % 16 evenly spaced pie-shaped quadrants
    
    % Store measures for each quadrant in nested cell array
    params.dwellQuad{j} = quad_measures.dwellQuad; % Dwell time
    params.pathLQuad{j} = quad_measures.pathLQuad; % distance travelled
    params.pathIVQuad{j} = quad_measures.pathIVQuad; % mean linear velocity
    params.numstopQuad{j} = quad_measures.numstopQuad; % number of stops
    params.stopQuad{j} = quad_measures.stopQuad; % stop duration 
    params.angVelQuad{j} = quad_measures.angVelQuad; % mean absolute angular velocity
    
    % Assess time near perimeter
    [params.dwellOutsideWall{j},params.centerDwell{j}] = OF.thigmotaxis(j,params,fr,.8);

 %% Build Occupancy Map
    [params.rateMap{j},map] = OF.occ_map(j,params,binsize,fr);
    
 %% Occupancy Map Measures (search area, home bases, home base area, all home base measures)
    params.searchArea{j} = OF.search_area(map);

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
    for hb = 1:size(params.HBBound{j},2)
        [ ~,~, tempC] = min_encl_ellipsoid( params.HBBound{j}{1, hb}(:,1),params.HBBound{j}{1, hb}(:,2));
        params.HBcenter{j}{1,hb} = [tempC(1,1),tempC(2,1)]; %x y coords
        
    end
    
    %Average Proximity of stops from hb center
    for hb = 1:size(params.HBcenter{j},2)
        tempHB = params.HBcenter{j}{1,hb};
        for i = 1:size(params.stops{j},2)
            
            temp = params.stops{j}{1,i};
            dist = sqrt((tempHB(1,1)-temp(:,1)).^2+(tempHB(1,2)-temp(:,2)).^2);
            disttemp(i,1) = nanmean(dist);
            
        end
        params.distHBstop{j}{1,hb} = nanmean(disttemp); %average distance from stops 
        params.closeHBstop{j}{1,hb} = sum(disttemp<25); %number of stops within 25cm from hb center
        params.distHBstopAll{j} = disttemp; clear disttemp tempHB dist temp
    end
    
    %Calculate distance measures between high occupancy coordinates centers
    temp = [];
    if size(params.HBcenter{j},2)>1
        for hb = 1:size(params.HBcenter{j},2)
            temp = [temp; params.HBcenter{j}{1,hb}];
        end
        params.avgHBdist{j} = nanmean(pdist(temp));
        params.maxHBdist{j} = max(pdist(temp));
        params.minHBdist{j} = min(pdist(temp));
    else
        params.avgHBdist{j} = NaN;
        params.maxHBdist{j} = NaN;
        params.minHBdist{j} = NaN;
    end
    clear temp
    
    %%calculating proximity of high occupancy coordinates center from the
    %%cue boundary
    if ~isempty(params.cueCoords{j})
        k = convhull(params.cueCM{j}(:,1),params.cueCM{j}(:,2));
        cueBoundary = [params.cueCM{j}(k,1),params.cueCM{j}(k,2)];
        
        for r = 1:size(params.HBcenter{j},2)
            distances = sqrt(sum(bsxfun(@minus, cueBoundary, [params.HBcenter{j}{1,r}(:,1),params.HBcenter{j}{1,r}(:,2)]).^2,2));
            params.HBdist2Cue{j}{1,r} = unique(distances(distances==min(distances))); %Find minimum distance from cue boundary to hb center
        end
    else
         k = convhull(params.cueCM{j-1}(:,1),params.cueCM{j-1}(:,2));
        cueBoundary = [params.cueCM{j-1}(k,1),params.cueCM{j-1}(k,2)];
        
        for r = 1:size(params.HBcenter{j},2)
            distances = sqrt(sum(bsxfun(@minus, cueBoundary, [params.HBcenter{j}{1,r}(:,1),params.HBcenter{j}{1,r}(:,2)]).^2,2));
            params.HBdist2Cue{j}{1,r} = unique(distances(distances==min(distances))); %Find minimum distance from cue boundary to hb center
        end
    end
    
%% Excursion Analysis 
   
    %This finds stops that occur in the home base boundary
    [startStop,endStop,~] = findgroups(params.stopIdx{j});
    
    for hb = 1:size(params.HBBound{j},2)
        endIn{hb} = inpolygon(back(endStop',1),back(endStop',2), params.HBBound{j}{1, hb}(:,1),params.HBBound{j}{1, hb}(:,2));
        startIn{hb} = inpolygon(back(startStop',1),back(startStop',2), params.HBBound{j}{1, hb}(:,1),params.HBBound{j}{1, hb}(:,2));
    end
    
    params.startInHB{j} = endIn; clear endIn
    params.endInHB{j} = startIn; clear startIn
    
 %find time when moving greater than or equal to 3cm/s for 1 sec 
    params.runIdx{j} = contiguousframes(params.pathIV{j}>=3,30);
    [startMotion,endMotion,~] = findgroups(params.runIdx{j});
    
    for ii = 1:length(startMotion)
        motion{ii} = [back(startMotion(ii):endMotion(ii),1),back(startMotion(ii):endMotion(ii),2)];
        timeMotion{ii} = size(motion{ii},1)/fr;
        tsMotion{ii} = ts(startMotion(ii):endMotion(ii),1);
    end
    
    params.runs{j} = motion;   clear motion
    params.timeMoving{j} = timeMotion; clear timeMotion
    params.tsMotion{j} = tsMotion; clear tsMotion
    params.NumRuns{j} = size(params.runs{j},2);
    
end
