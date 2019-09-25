% time_near_cue
% Ryan Harvey 2019; updated by LB 2019 to included entries/stop measures
load('d:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\params_V17.mat')
tic
for i=1:size(params.subID)
    
    if any(isnan(params.cueCM{i}))
        % rescale cue by 5 cm
        scaled_cuecords_x=rescale(params.cueCM{i-1}(:,1),min(params.cueCM{i-1}(:,1))-5,max(params.cueCM{i-1}(:,1))+5);
        scaled_cuecords_y=rescale(params.cueCM{i-1}(:,2),min(params.cueCM{i-1}(:,2))-5,max(params.cueCM{i-1}(:,2))+5);
    else
        % rescale cue by 5 cm
        scaled_cuecords_x=rescale(params.cueCM{i}(:,1),min(params.cueCM{i}(:,1))-5,max(params.cueCM{i}(:,1))+5);
        scaled_cuecords_y=rescale(params.cueCM{i}(:,2),min(params.cueCM{i}(:,2))-5,max(params.cueCM{i}(:,2))+5);
    end
    
    % convex hull to close shape of cue
    k = convhull(scaled_cuecords_x,scaled_cuecords_y);
    scaled_cuecords_x=scaled_cuecords_x(k);
    scaled_cuecords_y=scaled_cuecords_y(k);
    
    % Entries in Cue Zone
    
    tempIn=inpolygon(params.backCM{i,1}(:,1),params.backCM{i,1}(:,2),scaled_cuecords_x',scaled_cuecords_y');
    if sum(tempIn) == 0
        params.CueEntries(i)=0;
    else
        tempOut=contiguousframes(tempIn,60); %has to be inside of cue for at least 2 sec to count as entry
        [~,~,params.CueEntries(i)]=findgroups(tempOut);
    end
    
    % Retrieve Stops in Cue Zone
    
    stops = params.stops{i};
    
    for ii = 1:length(stops)
        
        temp_stop = stops{1,ii}(1,1:2);
        
        stop_in_cue(ii,1) = inpolygon(temp_stop(:,1),...
            temp_stop(:,2),scaled_cuecords_x',scaled_cuecords_y'); %Find number of times animals initiated a stop in the cue zone
        
        % compile timestamps of stop in cue
        if stop_in_cue (ii,1) == 1
            
            %sanity check
            plot(scaled_cuecords_x',scaled_cuecords_y'); hold on; plot(temp_stop(1,1),temp_stop(1,2),'*r')
            
            temp_ts = params.tsStop{i}{1,ii}(1,1); % get first ts for stop
            
            tsStops(ii,1) = temp_ts;
            
        else
            tsStops(ii,1) = NaN;
            
        end
        
    end
    
    params.CueStops(i)=nansum(stop_in_cue); %number of stops in cue zone
    
    % Bin stops that occur in cue zone over time
    
    if exist('tsStops','var')
        edges = [0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30];
        edges = edges*60; % convert edges into seconds
        
        params.bin_stops{i} = histcounts(tsStops,edges);
    else
        params.bin_stops{i} = zeros(1,16);
    end
    
    clear stop_in_cue tsStops
    
    % locate path in cue zone
    innose=inpolygon(params.noseCM{i}(:,1),params.noseCM{i}(:,2),scaled_cuecords_x,scaled_cuecords_y);
    inhead=inpolygon(params.headCM{i}(:,1),params.headCM{i}(:,2),scaled_cuecords_x,scaled_cuecords_y);
    
    % conjuctive nose or head in bound
    inbound=innose | inhead;
    
    hd=wrapTo360(fixNLXangle(rad2deg(atan2(params.noseCM{i}(:,2)-params.headCM{i}(:,2),...
        params.noseCM{i}(:,1)-params.headCM{i}(:,1))),round(0.1667*30)));
    
    bd=wrapTo360(fixNLXangle(rad2deg(atan2(params.headCM{i}(:,2)-params.backCM{i}(:,2),...
        params.headCM{i}(:,1)-params.backCM{i}(:,1))),round(0.1667*30)));
    
    
    bincenters=(0:6:360)+3;
    bincenters(end)=[];
    
    % hd in and out of zone
    params.overall_hd{i}=wrapTo360(hd');
    params.overall_hd_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(hd,0:6:360)',deg2rad(6));
    
    params.hd_out_zone{i}=wrapTo360(hd(~inbound)');
    params.hd_out_zone_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(hd(~inbound),0:6:360)',deg2rad(6));
    
    if isempty(hd(inbound))
        params.hd_in_zone{i}=NaN;
        params.hd_in_zone_mvl{i}=NaN;
    else
        params.hd_in_zone{i}=wrapTo360(hd(inbound)');
        params.hd_in_zone_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(hd(inbound),0:6:360)',deg2rad(6));
    end
    
    % body direction in and out of zone
    params.overall_body_dir{i}=wrapTo360(bd');
    params.overall_body_dir_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(bd,0:6:360)',deg2rad(6));
    
    params.body_dir_out_zone{i}=wrapTo360(bd(~inbound)');
    params.body_dir_out_zone_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(bd(~inbound),0:6:360)',deg2rad(6));
    
    if isempty(bd(inbound))
        params.body_dir_in_zone{i}=NaN;
        params.body_dir_in_zone_mvl{i}=NaN;
    else
        params.body_dir_in_zone{i}=wrapTo360(bd(inbound)');
        params.body_dir_in_zone_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(bd(inbound),0:6:360)',deg2rad(6));
    end
    
    % egocentric head bearing
    polyin = polyshape(scaled_cuecords_x,scaled_cuecords_y);
    [x,y] = centroid(polyin);
    
    head_bearing=wrapTo360(fixNLXangle(rad2deg(atan2(y-params.headCM{i}(:,2),...
        x-params.headCM{i}(:,1))),round(0.1667*30)));
    
    head_bearing=wrapTo360(head_bearing-hd);
    
    % egocentric bearing in and out of zone
    params.overall_head_bearing{i}=wrapTo360(head_bearing');
    params.overall_head_bearing_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(head_bearing,0:6:360)',deg2rad(6));
    
    params.head_bearing_out_zone{i}=wrapTo360(head_bearing(~inbound)');
    params.head_bearing_out_zone_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(head_bearing(~inbound),0:6:360)',deg2rad(6));
    
    if isempty(head_bearing(inbound))
        params.head_bearing_in_zone{i}=NaN;
        params.head_bearing_in_zone_mvl{i}=NaN;
    else
        params.head_bearing_in_zone{i}=wrapTo360(head_bearing(inbound)');
        params.head_bearing_in_zone_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(head_bearing(inbound),0:6:360)',deg2rad(6));
    end
    
    
    % egocentric body bearing
    polyin = polyshape(scaled_cuecords_x,scaled_cuecords_y);
    [x,y] = centroid(polyin);
    
    body_bearing=wrapTo360(fixNLXangle(rad2deg(atan2(y-params.backCM{i}(:,2),...
        x-params.backCM{i}(:,1))),round(0.1667*30)));
    
    body_bearing=wrapTo360(body_bearing-bd);
    
    % egocentric bearing in and out of zone
    params.overall_body_bearing{i}=wrapTo360(body_bearing');
    params.overall_body_bearing_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(body_bearing,0:6:360)',deg2rad(6));
    
    params.body_bearing_out_zone{i}=wrapTo360(body_bearing(~inbound)');
    params.body_bearing_out_zone_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(body_bearing(~inbound),0:6:360)',deg2rad(6));
    
    if isempty(body_bearing(inbound))
        params.body_bearing_in_zone{i}=NaN;
        params.body_bearing_in_zone_mvl{i}=NaN;
    else
        params.body_bearing_in_zone{i}=wrapTo360(body_bearing(inbound)');
        params.body_bearing_in_zone_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(body_bearing(inbound),0:6:360)',deg2rad(6));
    end
    
    % time in and out of zone
    params.time_in_zone{i}=sum(inbound)/30;
    params.time_out_zone{i}=sum(~inbound)/30;
    
end
toc