
% TESTING_HORSESHOE_LAP_PLOT

% root=data.HPCatn02.S20180715143806.root  ;
ncells=length(data.(root{1}).(root{2}).Spikes);
nsessions=size(data.(root{1}).(root{2}).events,2);
if any(data.(root{1}).(root{2}).linear_track_manage==1)
    nsessions=size(data.(root{1}).(root{2}).events,2)+1;
end

horseshoe=data.(root{1}).(root{2}).linear_track.nonlinearFrames;
for i=1:ncells
    figure;
    for ns=1:nsessions
        
        for lap=1:length(data.(root{1}).(root{2}).linear_track.right{1,i}.laps)
            horselapsR{lap,1}=horseshoe(ismember(horseshoe(:,1),...
                data.(root{1}).(root{2}).linear_track.right{1,i}.laps{lap}(:,1)),:);
        end
        
        for lap=1:length(data.(root{1}).(root{2}).linear_track.left{1,i}.laps)
            horselapsL{lap,1}=horseshoe(ismember(horseshoe(:,1),...
                data.(root{1}).(root{2}).linear_track.left{1,i}.laps{lap}(:,1)),:);
        end
    end
    spkframes=data.(root{1}).(root{2}).linear_track.right{1,i}.dataspks;
    spkts=spkframes(spkframes(:,6)==1,1);
    subplot(1,2,1)
    horselaps_togetherR=vertcat(horselapsR{:});
    plot(horselaps_togetherR(:,2),horselaps_togetherR(:,3),'.k');hold on
    [ts,idx]=unique(horselaps_togetherR(:,1));
    scatter(interp1(ts,horselaps_togetherR(idx,2),spkts),interp1(ts,horselaps_togetherR(idx,3),spkts),'Filled','r')
    
    spkframes=data.(root{1}).(root{2}).linear_track.left{1,i}.dataspks;
    spkts=spkframes(spkframes(:,6)==1,1);
    subplot(1,2,2)
    horselaps_togetherL=vertcat(horselapsL{:});
    plot(horselaps_togetherL(:,2),horselaps_togetherL(:,3),'.k');hold on
    [ts,idx]=unique(horselaps_togetherL(:,1));
    scatter(interp1(ts,horselaps_togetherL(idx,2),spkts),interp1(ts,horselaps_togetherL(idx,3),spkts),'Filled','r')
    
    clear horselapsL horselapsR
end
% I NEED TO PROPERLY INTERPOLATE SPIKE TIMES BETWEEN HORSESHOE, JUST LIKE
% IN createframes_w_spikebinary

% horselaps_togetherR=vertcat(horselapsR{:});
% figure;
% plot(horselaps_togetherR(:,2),horselaps_togetherR(:,3),'.k');hold on
% 
% [ts,idx]=unique(horselaps_togetherR(:,1));
% scatter(interp1(ts,horselaps_togetherR(idx,2),spkts),interp1(ts,horselaps_togetherR(idx,3),spkts),'Filled','r')
% 
% 
% horselaps_togetherL=vertcat(horselapsL{:});
% figure;
% plot(horselaps_togetherL(:,2),horselaps_togetherL(:,3),'.r');
% [ts,idx]=unique(horselaps_togetherL(:,1));
% scatter(interp1(ts,horselaps_togetherL(idx,2),spkts),interp1(ts,horselaps_togetherL(idx,3),spkts),'Filled','r')
% 

