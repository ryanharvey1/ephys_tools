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
fr = 30; %frame rate
binsize = 3; %to use for HB detection
upFactor = 15; %used in segmentimage - this upscales the heatmap to get a smooth boundary f

param_idx = params.subID; %serves as index for getParam
%Loop through subjects

for j = 1:size(param_idx,1)
    % Set maze parameters
    diameter = params.dia{j};
    
    % Establish position coordinates
    x = params.backCM{j}(:,1);
    y = params.backCM{j}(:,2);
    
    %Grab variables for easier use
    head = params.headCM{j};
    ts  = params.ts{j};
    
    % Calculate linear velocity, acceleration, and distance vector
    [params.velocity{j}, params.acceleration{j}, params.distance_vector{j}] = linear_motion(x,y,fr,3);
    
    
    %% Basic locomotor Measures (path length, velocity, number of stops/stop duration, time spent in outer/inner part of maze)
    
    params.path_length{j} = OF.path_length(x,y,vel,fr);
    
    
    % Compute Stopping behaviors
    stop_measures = OF.stops(x,y,ts,params.velocity{j},fr,epoch);
    
    params.stopIdx{j} = stop_measures.stopIdx;
    params.stops{j} = stop_measures.stops;
    params.timeStopped{j} = stop_measures.timeStopped;
    params.tsStop{j} = stop_measures.tsStop;
    params.NumStops{j} = stop_measures.NumStops;
    params.stopCenter{j} = stop_measures.stopCenter;
    
    % Assess time near perimeter
    [params.thigmotaxis{j}] = OF.thigmotaxis(x,y,fr,diameter,.8);
    
    %% Occupancy Map Measures (search area, home bases, home base area, all home base measures)
    [params.occupancy_map{j},map] = OF.occ_map(j,params,binsize,fr);
    
    % Search area
    params.searchArea{j} = OF.search_area(map);
    
    % High occupancy coordinates
    [~,~,home_base_x,home_base_y,params.fieldarea{j},params.upHBmap{j}] = segmentImage('map',params.rateMap{j},'upscalefac',upFactor); %Store area
    
    % Get home base metrics for each identified home base
    for f = 1:length(home_base_x)
        
        % rescale the home base coordinates
        [params.HBBound{j}(1,f),params.HBcoords{j}(1,f),params.HBcenter{j}(1,f),out_home,x_home,y_home] = ...
            OF.rescale_home_base(home_base_x{f},home_base_y{f},params.upHBmap{j},upFactor,diameter,fr);
        
        metrics = OF.home_base_metics(out_home,x,y,velocity,x_home,y_home,stopIdx,params.tsStop{j});
        
        % For homebases where animal didn't move majorily slow, pass on
        % saving metrics
        if metrics.HBclass > .8
            params.entries{j}{1,f} = NaN;
            params.hbOcc{j}(1,f) = NaN;
            params.hbVel{j}(1,f) = NaN;
            params.slowInHB{j}{1,f} = NaN;
            params.HBclass{j}{1,f} = NaN;
            params.HBstops{j}{1,f} = NaN;
            params.time2HB{j}{1,f} = NaN;
            params.HB_stop_dist{j}(1,f) = NaN;
            params.HB_close_stop{j}(1,f) = NaN;
            params.HB_stop_dist_vector{j}(1,f) = NaN;
            
        else
            
            % Average Proximity of stops from hb center
            [params.HB_stop_dist{j}(1,f),params.HB_close_stop{j}(1,f),params.HB_stop_dist_vector{j}(1,f)] = OF.stops_to_homebase(params.HBcenter{j}(1,f),params.stops{j});
            
            % unpack and save metrics
            params.entries{j}{1,f} = metrics.entries;
            params.hbOcc{j}(1,f) = metrics.hbOcc;
            params.hbVel{j}(1,f) = metrics.hbVel;
            params.slowInHB{j}{1,f} = metrics.slowInHB;
            params.HBclass{j}{1,f} = metrics.HBclass;
            params.HBstops{j}{1,f} = metrics.HBstops;
            params.time2HB{j}{1,f} = metrics.time2HB;
        end
        
    end
    
    
    % Calculate distance measures between high occupancy coordinates centers
        [params.HB_avg_dist{j},params.HB_max_dist{j},params.HB_min_dist{j}] = ...
            OF.distance_between_homebase(params.HBcenter{j});

    
        % calculating proximity of high occupancy coordinates center from the
        % cue boundary
        if ~isempty(params.cueCoords{j})
            params.HBdist2Cue{j} =  homebase_dist_from_cue(params.cueCM{j},params.HBcenter{j});
        else
            params.HBdist2Cue{j} =  homebase_dist_from_cue(params.cueCM{j - 1},params.HBcenter{j});
        end
    
    %% Excursion Analysis
    
    %This finds stops that occur in the home base boundary
    [startStop,endStop,~] = findgroups(params.stopIdx{j});
    
    for hb = 1:size(params.HBBound{j},2)
        endIn{hb} = inpolygon(x(endStop'),y(endStop'), params.HBBound{j}{1, hb}(:,1),params.HBBound{j}{1, hb}(:,2));
        startIn{hb} = inpolygon(x(startStop'),y(startStop'), params.HBBound{j}{1, hb}(:,1),params.HBBound{j}{1, hb}(:,2));
    end
    
    params.startInHB{j} = endIn; clear endIn
    params.endInHB{j} = startIn; clear startIn
    
    %find time when moving greater than or equal to 3cm/s for 1 sec
    params.runIdx{j} = contiguousframes(params.velocity{j}> = 3,30);
    [startMotion,endMotion,~] = findgroups(params.runIdx{j});
    
    for ii = 1:length(startMotion)
        motion{ii} = [x(startMotion(ii):endMotion(ii)),y(startMotion(ii):endMotion(ii))];
        timeMotion{ii} = size(motion{ii},1)/fr;
        tsMotion{ii} = ts(startMotion(ii):endMotion(ii),1);
    end
    
    params.runs{j} = motion;   clear motion
    params.timeMoving{j} = timeMotion; clear timeMotion
    params.tsMotion{j} = tsMotion; clear tsMotion
    params.NumRuns{j} = size(params.runs{j},2);
    
end
