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
load('params_V11.mat') %%loads table containing the raw coords and xMax/Min and yMax/Min
clearvars -except params

%% _________________________Compute Measures_______________________________



param_idx=params.PathName; %serves as index for getParam

%Loop through subjects
for j=1:size(param_idx,1)
    
    if ~contains(param_idx{j},'lgOF') % to bypass small open field for now LB March 2019
        continue
    end
    
    close all
    %% Maze Parameters;
    
    
    % compute total distance (pathL) and instantaneous velocity (pathIV)
    [ params.pathL{j},params.pathIV{j}] = compute_OFpathCalc(params.transcoords{j}(:,1),params.transcoords{j}(:,2),fr);
    
    
    %calculate verticies for quadrants
    quadrants=createZones([0,0],params.dia(j),'numQuad',16,'fig',0); %16 pie shaped segments
    
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
        [ segPL{ii},segPathIV{ii}] = compute_OFpathCalc(params.transcoords{j}(start(ii):ends(ii),1),...
            params.transcoords{j}(start(ii):ends(ii),2),fr);
        %         plot(params.transcoords{j}(start(ii):ends(ii),1),params.transcoords{j}(start(ii):ends(ii),2));hold on
    end
    params.runs{j}=motion;   clear motion
    params.segPL{j}=segPL; clear segPL
    params.segPathIV{j}=segPathIV; clear segPathIV
    
    %Compute circuity between segments
    for iii=1:size(params.runs{j},2)
        [params.circuity{j}(1,iii)]=circuity(params.runs{j}{iii}(:,1),params.runs{j}{iii}(:,2));
    end
    
    params.meanCirc{j}=nanmean(params.circuity{j}(1,:)); %compute average circuity
    
    
    params.NumRuns{j}= size(params.runs{j},2);
    
    %Creates bin edges for heatmap
    xedge=linspace(-(params.dia(j)/2),(params.dia(j)/2),round(params.dia(j)/binsize));
    yedge=linspace(-(params.dia(j)/2),(params.dia(j)/2),round(params.dia(j)/binsize));
    
    %Bin coordinates for heat map and apply guassian filter to smooth
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
    
    %Find number of stops with 25cm diameter of hb center
    params.stopsNearHB{j}=sum(params.proximity{j}<25);
    
    %Get path segments of inward bounds to HB
    
    
    %Get path segments of outward bounds from HB
end
