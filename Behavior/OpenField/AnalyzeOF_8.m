
%%
%%%%%%%%%%% Plot Paths with home base outline and cue location

load('params_V17.mat')
A = table2array(params);
param_idx=params.subID;

%% Get Primary home base measures
numHb=cellfun('size',params.hbOcc,2); 

for i=1:size(params.HBcenter,1)
    [~,occIdx(i,1)]=max(params.hbOcc{i});
    primary_hbOcc(i,1)=params.hbOcc{i}(1,occIdx(i,1));
    primary_hbEntry(i,1)=params.entries{i}{1,occIdx(i,1)};
    primary_hbDist2Cue(i,1)=params.HBdist2Cue{i}{1,occIdx(i,1)};
    primary_hbStops(i,1)=params.HBstops{i}{1,occIdx(i,1)};
    primary_hbvel(i,1)=params.hbVel{i}(1,occIdx(i,1));
    primary_area(i,1)=params.fieldarea{i}(1,occIdx(i,1));
    primary_distHBstop(i,1)=params.distHBstop{i}{1,occIdx(i,1)};
    primary_numCloseHBstops(i,1)=params.closeHBstop{i}{1,occIdx(i,1)};
    primary_hbDistBinary(i,1)=params.HBdist2Cue{i}{1,occIdx(i,1)}<=25;
    primary_time2HB(i,1)=params.time2HB{i}(1,occIdx(i,1));
end


%Create contingency Tables 
tg1=primary_hbDistBinary(contains(param_idx,'Tg') & contains(param_idx,'D1'),1);
wt1=primary_hbDistBinary(contains(param_idx,'WT') & contains(param_idx,'D1'),1);
tg2=primary_hbDistBinary(contains(param_idx,'Tg') & contains(param_idx,'D2'),1);
wt2=primary_hbDistBinary(contains(param_idx,'WT') & contains(param_idx,'D2'),1);

%Day 1 close/far
%row: tg,wt, col:close,far
day1=[sum(tg1) sum(tg1==0); sum(wt1) sum(wt1==0)];
[d1h,d1p,d1stats] = fishertest(day1);
%Day 1 close/far
%row: tg,wt, col:close,far
day2=[sum(tg2) sum(tg2==0); sum(wt2) sum(wt2==0)];
[d2h,d2p,d2stats] = fishertest(day2);

timeInCue=cell2mat(params.time_in_zone);

vars={'Subject','day','group','NumHBs','Occupancy','Entries',...
    'Distance2Cue','NumStops','AvgVelocity','area','avgStopDist',...
    'numCloseStops','clos2cue','time2hb','timeInCue'}; %varnames act as headers for hbMeas2Python
day=repmat({'D1';'D2'},size(params.HBcenter,1)/2,1); %create day variable 
group=[repmat({'Tg'},size(params.HBcenter,1)/2,1);...%create grouping variable
    repmat({'Wt'},size(params.HBcenter,1)/2,1)];
id=split(params.subID,'_'); %create subID as a within-subjects factor
[~,~,subject]=unique(id(:,1)); %create subID as a within-subjects factor 

hbMeas2Python=table(subject,day,group,numHb,primary_hbOcc,...
    primary_hbEntry, primary_hbDist2Cue,...
    primary_hbStops, primary_hbvel,primary_area,...
    primary_distHBstop,primary_numCloseHBstops,...
    primary_hbDistBinary,primary_time2HB,timeInCue,'VariableNames',vars); %make table

writetable(hbMeas2Python,...
    'D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\hbData.csv'); %save data

%% Proximity of primary hb day 1 versus primary hb day 2
row=1;
for i=1:2:size(params,1)
    primaryHBdist(row,1)=sqrt((params.HBcenter{i,1}{1,occIdx(i,1)}(1,1)-params.HBcenter{i+1,1}{1,occIdx(i+1,1)}(1,1))^2+...
        (params.HBcenter{i,1}{1,occIdx(i,1)}(1,2)-params.HBcenter{i+1,1}{1,occIdx(i+1,1)}(1,2))^2);
    row=row+1;
end

%% Secondary home bases 

for i=1:size(params.HBcenter,1)
    [~,occIdx]=sort(params.hbOcc{i},'descend');
    if length(occIdx)>=2
    secondary_hbOcc(i,1)=params.hbOcc{i}(1,occIdx(2));
    secondary_hbEntry(i,1)=params.entries{i}{1,occIdx(2)};
    secondary_hbDist2Cue(i,1)=params.HBdist2Cue{i}{1,occIdx(2)};
    secondary_hbStops(i,1)=params.HBstops{i}{1,occIdx(2)};
    secondary_hbvel(i,1)=params.hbVel{i}(1,occIdx(2));
    secondary_area(i,1)=params.fieldarea{i}(1,occIdx(2));
    secondary_distHBstop(i,1)=params.distHBstop{i}{1,occIdx(2)};
    secondary_numCloseHBstops(i,1)=params.closeHBstop{i}{1,occIdx(2)};
    secondary_hbDistBinary(i,1)=params.HBdist2Cue{i}{1,occIdx(2)}<=25;
    secondary_time2HB(i,1)=params.time2HB{i}(1,occIdx(2));
    else
    secondary_hbOcc(i,1)=NaN;
    secondary_hbEntry(i,1)=NaN;
    secondary_hbDist2Cue(i,1)=NaN;
    secondary_hbStops(i,1)=NaN;
    secondary_hbvel(i,1)=NaN;
    secondary_area(i,1)=NaN;
    secondary_distHBstop(i,1)=NaN;
    secondary_numCloseHBstops(i,1)=NaN;
    secondary_hbDistBinary(i,1)=NaN;
    end
    
end

figure; 
subplot(1,2,1)
plotspread_wrapper(cell2mat(wt1),cell2mat(tg1),{'WT','Tg'})
title('Time to first stop in Primary HB - Cue Present')
ylim([0 max([cell2mat(wt1);cell2mat(tg1);cell2mat(wt2);cell2mat(tg2)])])
ylabel('time (s)')
subplot(1,2,2)
plotspread_wrapper(cell2mat(wt2),cell2mat(tg2),{'WT','Tg'})
title('Time to first stop in Primary HB - Cue Absent')
ylim([0 max([cell2mat(wt1);cell2mat(tg1);cell2mat(wt2);cell2mat(tg2)])])
ylabel('time (s)')


figure; 
subplot(1,2,1)
plotspread_wrapper(secondary_hbDist2Cue(25:2:end,1),secondary_hbDist2Cue(1:2:24,1),{'WT','Tg'})
title('Day 1 - HB dist from cue')
subplot(1,2,2)
plotspread_wrapper(secondary_hbDist2Cue(26:2:end,1),secondary_hbDist2Cue(2:2:24,1),{'WT','Tg'})
title('Day 2 - HB dist from cue')

fig=figure; 
fig.Color=[1 1 1];
scatter(primary_hbDist2Cue(1:24,1),secondary_hbDist2Cue(1:24,1),'r','filled')
hold on; 
scatter(primary_hbDist2Cue(25:end,1),secondary_hbDist2Cue(25:end,1),'k','filled')


%% Create contingency Tables 
tg1=secondary_hbDistBinary(contains(param_idx,'Tg') & contains(param_idx,'D1'),1);
wt1=secondary_hbDistBinary(contains(param_idx,'WT') & contains(param_idx,'D1'),1);
tg2=secondary_hbDistBinary(contains(param_idx,'Tg') & contains(param_idx,'D2'),1);
wt2=secondary_hbDistBinary(contains(param_idx,'WT') & contains(param_idx,'D2'),1);

%Day 1 close/far
%row: tg,wt, col:close,far
day1=[nansum(tg1) nansum(tg1==0); nansum(wt1) nansum(wt1==0)];
[d1h,d1p,d1stats] = fishertest(day1)
%Day 1 close/far
%row: tg,wt, col:close,far
day2=[nansum(tg2) nansum(tg2==0); nansum(wt2) nansum(wt2==0)];
[d2h,d2p,d2stats] = fishertest(day2)

%% Evaluate cue quadrant measures 

tg1=vertcat(params.dwellQuad{contains(param_idx,'Tg') & contains(param_idx,'D1')});
wt1=vertcat(params.dwellQuad{contains(param_idx,'WT') & contains(param_idx,'D1')});
tg2=vertcat(params.dwellQuad{contains(param_idx,'Tg') & contains(param_idx,'D2')});
wt2=vertcat(params.dwellQuad{contains(param_idx,'WT') & contains(param_idx,'D2')});

subplot(1,2,1)
plotspread_wrapper(sum(wt1(:,5:7),2),sum(tg1(:,5:7),2),{'WT','Tg'})
title('Cue quadrant dwell time')

subplot(1,2,2)
plotspread_wrapper(sum(wt2(:,5:7),2),sum(tg2(:,5:7),2),{'WT','Tg'})
title('Cue quadrant dwell time')

day=[ones(size([sum(tg1(:,5:7),2);sum(wt1(:,5:7),2)],1),1);repmat(2,size([sum(tg1(:,5:7),2);sum(wt1(:,5:7),2)],1),1)];
group=[ones(size(tg1,1),1); zeros(size(wt1,1),1); ones(size(tg2,1),1); zeros(size(wt2,1),1)];
dwellQuad_stats=[sum(tg1(:,5:7),2); sum(wt1(:,5:7),2); sum(tg2(:,5:7),2); sum(wt2(:,5:7),2)];

tg1=abs(vertcat(params.angVelQuad{contains(param_idx,'Tg') & contains(param_idx,'D1')}));
wt1=abs(vertcat(params.angVelQuad{contains(param_idx,'WT') & contains(param_idx,'D1')}));
tg2=abs(vertcat(params.angVelQuad{contains(param_idx,'Tg') & contains(param_idx,'D2')}));
wt2=abs(vertcat(params.angVelQuad{contains(param_idx,'WT') & contains(param_idx,'D2')}));


subplot(1,2,1)
plotspread_wrapper(sum(wt1(:,5:7),2),sum(tg1(:,5:7),2),{'WT','Tg'})
title('Ang Vel Cue Quadrant')

subplot(1,2,2)
plotspread_wrapper(sum(wt2(:,5:7),2),sum(tg2(:,5:7),2),{'WT','Tg'})
title('Ang Vel Cue Quadrant')

AngQuad_stats=[nanmean(tg1(:,5:7),2); nanmean(wt1(:,5:7),2); nanmean(tg2(:,5:7),2); nanmean(wt2(:,5:7),2)];

tg1=vertcat(params.pathIVQuad{contains(param_idx,'Tg') & contains(param_idx,'D1')});
wt1=vertcat(params.pathIVQuad{contains(param_idx,'WT') & contains(param_idx,'D1')});
tg2=vertcat(params.pathIVQuad{contains(param_idx,'Tg') & contains(param_idx,'D2')});
wt2=vertcat(params.pathIVQuad{contains(param_idx,'WT') & contains(param_idx,'D2')});


subplot(1,2,1)
plotspread_wrapper(sum(wt1(:,5:7),2),sum(tg1(:,5:7),2),{'WT','Tg'})
title('Linear Vel Cue Quadrant')

subplot(1,2,2)
plotspread_wrapper(sum(wt2(:,5:7),2),sum(tg2(:,5:7),2),{'WT','Tg'})
title('Linear Vel Cue Quadrant')

linearQuad_stats=[nanmean(tg1(:,5:7),2); nanmean(wt1(:,5:7),2); nanmean(tg2(:,5:7),2); nanmean(wt2(:,5:7),2)];

tg1=vertcat(params.numstopQuad{contains(param_idx,'Tg') & contains(param_idx,'D1')});
wt1=vertcat(params.numstopQuad{contains(param_idx,'WT') & contains(param_idx,'D1')});
tg2=vertcat(params.numstopQuad{contains(param_idx,'Tg') & contains(param_idx,'D2')});
wt2=vertcat(params.numstopQuad{contains(param_idx,'WT') & contains(param_idx,'D2')});


subplot(1,2,1)
plotspread_wrapper(sum(wt1(:,5:7),2),sum(tg1(:,5:7),2),{'WT','Tg'})
title('Num Stops Cue Quadrant')

subplot(1,2,2)
plotspread_wrapper(sum(wt2(:,5:7),2),sum(tg2(:,5:7),2),{'WT','Tg'})
title('Num Stops Cue Quadrant')

numStopQuad_stats=[sum(tg1(:,5:7),2); sum(wt1(:,5:7),2); sum(tg2(:,5:7),2); sum(wt2(:,5:7),2)];

tg1=vertcat(params.pathLQuad{contains(param_idx,'Tg') & contains(param_idx,'D1')});
wt1=vertcat(params.pathLQuad{contains(param_idx,'WT') & contains(param_idx,'D1')});
tg2=vertcat(params.pathLQuad{contains(param_idx,'Tg') & contains(param_idx,'D2')});
wt2=vertcat(params.pathLQuad{contains(param_idx,'WT') & contains(param_idx,'D2')});

pathLQuad_stats=[sum(tg1(:,5:7),2); sum(wt1(:,5:7),2); sum(tg2(:,5:7),2); sum(wt2(:,5:7),2)];
subplot(1,2,1)
plotspread_wrapper(sum(wt1(:,5:7),2),sum(tg1(:,5:7),2),{'WT','Tg'})
title('Path Length Cue Quadrant')

subplot(1,2,2)
plotspread_wrapper(sum(wt2(:,5:7),2),sum(tg2(:,5:7),2),{'WT','Tg'})
title('Path LengthCue Quadrant')

vars={'Subject','day','group','dwell','angVel','linVel','numStops','pathL'};

sub=1:1:24;

subject=[sub';sub']; clear sub

quadMeas2Python=table(subject,day,group,dwellQuad_stats,AngQuad_stats,linearQuad_stats, numStopQuad_stats,...
    pathLQuad_stats,'VariableNames',vars);

writetable(quadMeas2Python,'D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\quadData.csv');
