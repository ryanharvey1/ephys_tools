working...
function decoding_preprocessing(data,event)
    linear_track=data.linear_track;
    spikesID=data.spikesID;
    
    load(uniquesessions{i,1},'spikesID','linear_track')
    
    ensemble_rows=find(contains(uniqueids(:,1),uniquesessions{i,1}));

    pos=[linear_track.right{1,1}.dataspks(linear_track.right{1,1}.dataspks(:,6)==0,:);...
        linear_track.left{1,1}.dataspks(linear_track.left{1,1}.dataspks(:,6)==0,:)];
    [~,I]=sort(pos(:,1));
    pos=pos(I,:);
    
    for k=1:length(ensemble_rows)
        cells=find(contains(spikesID.TetrodeNum,uniqueids{ensemble_rows(k),2})...
            & ismember(spikesID.CellNum,str2double(uniqueids(ensemble_rows(k),3))))';
        spike_times{k,1}=sort([linear_track.right{1,cells}.dataspks(linear_track.right{1,cells}.dataspks(:,6)==1,1);...
            linear_track.left{1,cells}.dataspks(linear_track.left{1,cells}.dataspks(:,6)==1,1)]);
    end
    pos_times=pos(:,1);
    vel=pos(:,5);
    head_angle=pos(:,4);
    pos(:,[1,4:end])=[];
    
    
    save(['D:\Projects\HPCatn\DecodingData\',uniquesessions{i,1}],'spike_times','pos','pos_times','vel','head_angle')
    clear spike_times pos pos_times vel head_angle

end



