% time_near_cue
% Ryan Harvey 2019
load('D:\Projects\TG_largeOpenField\params_V12.mat')

for i=1:size(params.subID)
    
    if any(isnan(params.cueCM{i}))
        continue
    end
    
    % rescale cue by 5 cm
    scaled_cuecords_x=rescale(params.cueCM{i}(:,1),min(params.cueCM{i}(:,1))-5,max(params.cueCM{i}(:,1))+5);
    scaled_cuecords_y=rescale(params.cueCM{i}(:,2),min(params.cueCM{i}(:,2))-5,max(params.cueCM{i}(:,2))+5);
    
    % convex hull to close shape of cue
    k = convhull(scaled_cuecords_x,scaled_cuecords_y);
    scaled_cuecords_x=scaled_cuecords_x(k);
    scaled_cuecords_y=scaled_cuecords_y(k);

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
    params.overall_hd{i}=wrapTo360(rad2deg(circ_mean(deg2rad(hd'))));
    params.overall_hd_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(hd,0:6:360)',deg2rad(6));

    params.hd_out_zone{i}=wrapTo360(rad2deg(circ_mean(deg2rad(hd(~inbound)'))));
    params.hd_out_zone_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(hd(~inbound),0:6:360)',deg2rad(6));

    if isempty(hd(inbound))
        params.hd_in_zone{i}=NaN;
        params.hd_in_zone_mvl{i}=NaN;
    else
        params.hd_in_zone{i}=wrapTo360(rad2deg(circ_mean(deg2rad(hd(inbound)'))));
        params.hd_in_zone_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(hd(inbound),0:6:360)',deg2rad(6));
    end
    
    % body direction in and out of zone
    params.overall_body_dir{i}=wrapTo360(rad2deg(circ_mean(deg2rad(bd'))));
    params.overall_body_dir_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(bd,0:6:360)',deg2rad(6));

    params.body_dir_out_zone{i}=wrapTo360(rad2deg(circ_mean(deg2rad(bd(~inbound)'))));
    params.body_dir_out_zone_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(bd(~inbound),0:6:360)',deg2rad(6));

    if isempty(bd(inbound))
        params.body_dir_in_zone{i}=NaN;
        params.body_dir_in_zone_mvl{i}=NaN;
    else
        params.body_dir_in_zone{i}=wrapTo360(rad2deg(circ_mean(deg2rad(bd(inbound)'))));
        params.body_dir_in_zone_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(bd(inbound),0:6:360)',deg2rad(6));
    end
    
    % egocentric head bearing 
    polyin = polyshape(scaled_cuecords_x,scaled_cuecords_y);
    [x,y] = centroid(polyin);
    
    head_bearing=wrapTo360(fixNLXangle(rad2deg(atan2(y-params.headCM{i}(:,2),...
        x-params.headCM{i}(:,1))),round(0.1667*30)));
    
    head_bearing=wrapTo360(head_bearing-hd);
    
    % egocentric bearing in and out of zone
    params.overall_head_bearing{i}=wrapTo360(rad2deg(circ_mean(deg2rad(head_bearing'))));
    params.overall_head_bearing_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(head_bearing,0:6:360)',deg2rad(6));

    params.head_bearing_out_zone{i}=wrapTo360(rad2deg(circ_mean(deg2rad(head_bearing(~inbound)'))));
    params.head_bearing_out_zone_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(head_bearing(~inbound),0:6:360)',deg2rad(6));

    if isempty(head_bearing(inbound))
        params.head_bearing_in_zone{i}=NaN;
        params.head_bearing_in_zone_mvl{i}=NaN;
    else
        params.head_bearing_in_zone{i}=wrapTo360(rad2deg(circ_mean(deg2rad(head_bearing(inbound)'))));
        params.head_bearing_in_zone_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(head_bearing(inbound),0:6:360)',deg2rad(6));
    end
    
    
    % egocentric body bearing
    polyin = polyshape(scaled_cuecords_x,scaled_cuecords_y);
    [x,y] = centroid(polyin);
    
    body_bearing=wrapTo360(fixNLXangle(rad2deg(atan2(y-params.backCM{i}(:,2),...
        x-params.backCM{i}(:,1))),round(0.1667*30)));
    
    body_bearing=wrapTo360(body_bearing-bd);
    
    % egocentric bearing in and out of zone
    params.overall_body_bearing{i}=wrapTo360(rad2deg(circ_mean(deg2rad(body_bearing'))));
    params.overall_body_bearing_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(body_bearing,0:6:360)',deg2rad(6));

    params.body_bearing_out_zone{i}=wrapTo360(rad2deg(circ_mean(deg2rad(body_bearing(~inbound)'))));
    params.body_bearing_out_zone_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(body_bearing(~inbound),0:6:360)',deg2rad(6));

    if isempty(body_bearing(inbound))
        params.body_bearing_in_zone{i}=NaN;
        params.body_bearing_in_zone_mvl{i}=NaN;
    else
        params.body_bearing_in_zone{i}=wrapTo360(rad2deg(circ_mean(deg2rad(body_bearing(inbound)'))));
        params.body_bearing_in_zone_mvl{i}=circ_r(deg2rad(bincenters)',histcounts(body_bearing(inbound),0:6:360)',deg2rad(6));
    end
    
    % time in and out of zone
    params.time_in_zone{i}=sum(inbound)/30;
    params.time_out_zone{i}=sum(~inbound)/30;

end