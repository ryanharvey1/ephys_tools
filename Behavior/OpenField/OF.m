classdef OF
    %OF consists of functions used to assess locomotion and home base
    %   behavior during exploration of an open field.
    
    %   Functions were created to work with 'params' a subject x measure table created from
    %   OF_preprocess.m in ephys_tools > Behavior > OpenField.
    
    properties
        Property1
    end
    
    methods(Static)
        
        function [length,velocity,mean_vel,acceleration,mean_accel,distance_vector] = path_measures(j,params,fr)
            % Builds an occupancy map for maze. Assums params table.
            % inputs:
            %   j: subject index
            %   params: table made from OF_preprocess
            %   fr: frame rate
            % output:
            %   length: total path length in cm for active motion
            %   velocity: vector of instantaneous velocity (n - 1 length of
            %       path vector)
            %   acceleration: vector of instantanous acceleration (n - 1
            %       length of path vector).
            %   distance_vector: vector of distance between coordinates  (n - 1
            %       length of path vector).
            %   mean_vel: average velocity during motion.
            %   mean_accel: mean absolute acceleration when animal is
            %       moving over 3cm/s.
            
            back = params.backCM{j};
            x = back(:,1); y = back(:,2);
            
            %squared distance between consecutive points
            sqrXDiff = (diff(x)).^2;
            sqrYDiff = (diff(y)).^2;
            
            %distance formula
            distance_vector = sqrt(sqrXDiff + sqrYDiff);
            velocity = (distance_vector)*fr; %instanteous velocity
            acceleration = gradient(velocity,1/fr); %instanteous acceleration
            
            velocity = smoothdata(velocity,'movmedian',fr*.8); %smoothed with moving median window over less than 1 second
            
            % Summary Path Measures
            length = sum(distance_vector(velocity>=3,1));%Total Path length for points greater than 3cm/s
            
            % Mean Velocity during movement
            mean_vel = nanmean(velocity(velocity>=3,1));
            
            % Mean absolute acceleration
            accel_idx = acceleration ~= 0 ;
            run_idx = velocity >= 3;
            mean_accel = nanmean(abs(acceleration(run_idx & accel_idx,1)));
        end
        
        function stop_measures = stops(j,params,fr)
            % Finds when animal stops ( < 3cm/s velocity for at least 1
            % second) and computes stop features.
            % inputs:
            %       j: subject index
            %       params: table made from OF_preprocess
            %       fr: frame rate
            % outputs:
            %       stop_measures: structure containing outcome measures.
            %           stopIdx: logical index of points where rat is
            %                    stopped.
            %           stops: cell array of stop xy coordinates
            %           timeStopped: time (seconds) of a stop
            %           tsStop: cell array of time stamps corresponding to stop
            %           NumStops: number of stops made
            %
            back = params.backCM{j};
            ts=params.ts{j};
            x = back(:,1); y = back(:,2);
            %squared distance between consecutive points
            sqrXDiff = (diff(x)).^2;
            sqrYDiff = (diff(y)).^2;
            
            %distance formula
            distance_vector = sqrt(sqrXDiff + sqrYDiff);
            velocity = (distance_vector)*fr; %instanteous velocity
            velocity = smoothdata(velocity,'movmedian',fr*.8); %smoothed with moving median window over less than 1 second
            
            stop_measures.stopIdx=contiguousframes(velocity < 3,30);
            [startStop,endStop,~]=findgroups(stop_measures.stopIdx);
            
            %     This finds coords for stopping
            for ii=1:length(startStop)
                motionless{ii}=[back(startStop(ii):endStop(ii),1),back(startStop(ii):endStop(ii),2)];
                timeMotionless{ii}=size(motionless{ii},1)/fr;
                tsStop{ii}=ts(startStop(ii):endStop(ii),1);
            end
            
            stop_measures.stops=motionless;
            stop_measures.timeStopped=timeMotionless;
            stop_measures.tsStop=tsStop;
            
            %     Create Number of Stops
            stop_measures.NumStops = size(stop_measures.stops,2);
            
            %Find center of mass for stops.
            for ii=1:size(stop_measures.stops,2)
                temp=[stop_measures.stops{1,ii}(:,1),stop_measures.stops{1,ii}(:,2)];
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
            
            stop_measures.stopCenter=stop;
            
        end
        
        function quad_measures = quadrant(j,params,fr,numQuad)
            % Finds when animal stops ( < 3cm/s velocity for at least 1
            % second) and computes stop features.
            % inputs:
            %       j: subject index
            %       params: table made from OF_preprocess
            %       fr: frame rate
            % outputs:
            %       quad_measures: structure containing outcome measures.
            %           stopIdx: logical index of points where rat is
            %                    stopped.
            %           stops: cell array of stop xy coordinates
            %           timeStopped: time (seconds) of a stop
            %           tsStop: cell array of time stamps corresponding to stop
            %           NumStops: number of stops made
            %
            
            back = params.backCM{j};
            head = params.headCM{j};
            nose = params.noseCM{j};
            
            ts = params.ts{j};
            x = back(:,1); y = back(:,2);
            
            %squared distance between consecutive points
            sqrXDiff = (diff(x)).^2;
            sqrYDiff = (diff(y)).^2;
            
            %distance formula
            distance_vector = sqrt(sqrXDiff + sqrYDiff);
            velocity = (distance_vector)*fr; %instanteous velocity
            velocity = smoothdata(velocity,'movmedian',fr*.8); %smoothed with moving median window over less than 1 second
            
            clear sqrXDiff sqrYDiff
            
            %Stops
            stopIdx=contiguousframes(velocity < 3,30);
            [startStop,endStop,~]=findgroups(stopIdx);
            
            %     This finds coords for stopping
            for ii=1:length(startStop)
                motionless{ii}=[back(startStop(ii):endStop(ii),1),back(startStop(ii):endStop(ii),2)];
            end
            
            stops=motionless;
            
            clear startStop endStop
            
            %Find center of mass for stops.
            for ii=1:size(stops,2)
                temp=[stops{1,ii}(:,1),stops{1,ii}(:,2)];
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
            
            stopCenter=stop;
            
            clear stop
            
            HD=wrapTo360(fixNLXangle(rad2deg(atan2(head(:,2)-nose(:,2),...
                head(:,1)-nose(:,1))),round(0.1667*30)));
            angVel=insta_angvel(HD',fr);
            
            
            %calculate verticies for quadrants
            quadrants=createZones([0,0],params.dia{j},'numQuad',numQuad,'fig',0); %16 pie shaped segments
            
            
            %Calculate dwell time per quadrant
            for i=1:size(quadrants,1)-1
                tempXin = [0 quadrants(i,1) quadrants(i+1,1)];
                tempYin = [0 quadrants(i,2) quadrants(i+1,2)];
                [in,~] = inpolygon(back(:,1),back(:,2),tempXin,tempYin);
                quad_measures.dwellQuad{j}(1,i) = sum(in)/fr; %
                quad_measures.pathLQuad{j}(1,i) = sum(distance_vector(in(1:end-1) & velocity>3 ,1)); %mean distance when animal is in quadrant and moving >3cm/s
                quad_measures.pathIVQuad{j}(1,i) = nanmean(velocity(in(1:end-1),1)); %mean linear velocity
                quad_measures.numstopQuad{j}(1,i) = sum(inpolygon(stopCenter{j}(:,1),stopCenter{j}(:,2),tempXin,tempYin));
                quad_measures.stopQuad{j}(1,i) = sum(in(1:end-1,:) & stopIdx{j})/fr; %stop duration
                quad_measures.angVelQuad{j}(1,i) = nanmean(abs(angVel(in(1:end-1,:))));
                clear tempXin tempYin in
            end
            
        end
        
        function [occ,map] = occ_map(j,params,binsize,fr)
            % Builds an occupancy map for maze. Assums params table.
            % inputs:
            %   j: subject index
            %   params: table made from OF_preprocess
            %   binsize: how big each bin should be in cm
            %   fr: frame rate
            % output:
            %   occ: smoothed occupancy map
            %   map: raw occupancy map
            
            %Creates bin edges for heatmap
            xedge=linspace(-(params.dia{j}/2),(params.dia{j}/2),round(params.dia{j}/binsize));
            yedge=linspace(-(params.dia{j}/2),(params.dia{j}/2),round(params.dia{j}/binsize));
            
            %Bin coordinates for heat map and apply guassian filter to smooth
            back = params.backCM{j};
            [map] = histcounts2(back(:,1),back(:,2),xedge,yedge);
            map=flipud(imrotate(map,90));
            map=map/fr;
            occ = imgaussfilt(map, 1.5);
        end
        
        function [out,in] = thigmotaxis(j,params,fr,center_proportion)
            % Builds an occupancy map for maze. Assums params table.
            % inputs:
            %   j: subject index
            %   params: table made from OF_preprocess
            %   center_proportion: The inner proportion of the maze (.8
            %       would allow for thigmotaxis boundary of outer 20% of maze.
            %   fr: frame rate
            % output:
            %   occ: smoothed occupancy map

            back = params.backCM{j};
            
            %Create annuli for center and outter edge of maze
            outsideDwell=createZones([0,0],params.dia{j},'type','annulus','fig',0,'annulusSize',center_proportion); %default annulus size is 80% for createZones
            centerDwell=createZones([0,0],params.dia{j},'type','annulus','fig',0,'annulusSize',center_proportion); %default annulus size is 80% for createZones
            
            %Calculate dwell time for outter edge
            [in,~]=inpolygon(back(:,1),back(:,2),outsideDwell(:,1),outsideDwell(:,2));
            out=sum(~in)/fr;
            
            clear in outsideDwell
            
            %Calculate dwell time for center maze
            [in,~]=inpolygon(back(:,1),back(:,2),centerDwell(:,1),centerDwell(:,2));
            in=sum(in)/fr;

        end
        
        function sa = search_area(map)
            % computes proportion of maze area occupied by animal given
            % occupancy map.
            % inputs:
            %   map: raw occpancy map (nxn grid of binned occupancy
            %   weighted by frame rate)
            
            % output:
            %   sa: search area as a numeric proportion (e.g. .10 = 10%
            %   of maze searched).
            
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
            
            sa=sum(sum(map>0))/sum(sum(~isnan(map))); %Calculate the proportion of bins occupied by animal.
            
        end
        
    end
end

