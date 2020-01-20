
%% This script organizing data for the TgF344-AD open field study conducted by Berkowitz et al. 
%  Most data are exported as csv files for further analysis in Python.
%  However, some measures are directly analyzed below. 

load('D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\params_V19');
param_idx=params.subID;

%% Compile single numeric Measures for whole trial analysis (#1)

% These measures correspond to behaviors summarized across the entire trial
% pathL: total path length when animal is moving. 
% searchArea: proportion of maze the animal visited. 
% runSpeed: average linear running speed. 
% numStops: number of instances when the animal decreased their speed to
%           <3cm/s for at least 1 second. 
% CueEntries: number of entries in cue zone during probe
% CueStops: number of stops in cue zone during probe
%
% Measures are compiled to a table and exported as csv for analysis in
% Python. 

%Here we'll define the variable names. 
vars={'Subject','day','group','pathL','searchArea','runSpeed','runAngVel','stopAngVel'...
    ,'numRuns','numStops','outside_hb_stop','CueEntries','CueStops','runAcell'};

%Lets pull the variables from the table. 
pathL = cell2mat(params.pathL); 
searchArea = cell2mat(params.searchArea);
runSpeed = cell2mat(params.meanLinVel_run);
runAcell = params.med_lin_accel; 
runAngVel = cell2mat(params.meanAngVel_run);
stopAngVel = cell2mat(params.meanAngVel_stop);
numRuns = cell2mat(params.NumRuns);
numStops = cell2mat(params.NumStops); 
outside_hb_stop = cell2mat(params.outside_hb_stops);
CueEntries = params.CueEntries;
CueStops = params.CueStops;


%Now we'll create variables for levels of independent variables as well as
%make the unique subject id. 
day=repmat({'D1';'D2'},size(params.time_in_zone,1)/2,1); %create day variable 
group=[repmat({'Tg'},size(params.time_in_zone,1)/2,1);...%create grouping variable
    repmat({'Wt'},size(params.time_in_zone,1)/2,1)];
id=split(params.subID,'_'); %create subID as a within-subjects factor
[~,~,subject]=unique(id(:,1)); %create subID as a within-subjects factor 

%We can now create a table with only our new vars 
wholeTrial2Python=table(subject,day,group,pathL,searchArea,...
    runSpeed,runAngVel,stopAngVel,numRuns,numStops,outside_hb_stop,...
    CueEntries,CueStops,runAcell,'VariableNames',vars); %make table

writetable(wholeTrial2Python,...
    'D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\wholeTrial_measures.csv'); %save data

%% Distance of Runs and Circuity of Runs (Fig 1 E, D)

for i = 1:size(param_idx,1)

    median_circuity(i,1) = nanmedian(cell2mat(params.seg_circuity{i}));
    median_timeMove(i,1) = nanmedian(cell2mat(params.timeMoving{i}));
    median_timeStop(i,1) = nanmedian(cell2mat(params.timeStopped{i}));
    median_segDist(i,1) = nanmedian(cell2mat(params.seg_distance{i}));
    median_seg_stop_out(i,1) = nanmedian(params.stop_out_hb_dur{i});

end

vars = {'subject','day','group','seg_circ','seg_dur','stop_dur','seg_dist','seg_dur_out'};

seg2python = table(subject,day,group,median_circuity,...
    median_timeMove,median_timeStop,median_segDist,...
    median_seg_stop_out,'VariableNames',vars);

writetable(seg2python,...
    'D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\segment_measures.csv'); %save data

% eCDF of Circuity for all animals Fig 2 

param_idx=params.subID;
for i=1:2:48
    [f,x] = ecdf(cell2mat(params.seg_circuity{i})');
    if contains(param_idx(i),'Tg')
        plot(x,f,'LineWidth',3,'Color','r');
    else
        plot(x,f,'LineWidth',3,'Color',[.2 .2 .2]);
    end
    hold on;
end
grid on

ylabel('Proportion')
xlabel('Circuity')


%% Get Primary home base measures lb 9/18/19 Add mean stop duration in home base (#3) 


%Get logical of whether moving slow critera was met
for i=1:size(params.HBcenter,1)
    for ii = 1:size(params.hbOcc{i},2)
        move_slow{i,1}(1,ii) = cell2mat(params.HBclass{i}(1,ii)) > .75;
    end
end

%num of HB that meet above 
for i = 1 : length (params.hbOcc)
    numHB(i,1)= sum(move_slow{i});
end

for i=1:size(params.HBcenter,1)
    [~,occIdx(i,1)]=max(params.hbOcc{i}(move_slow{i})); %max given moving slow most of time;
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
    primary_hbStopDuration(i,1) = cell2mat(params.slowInHB{i}(1,occIdx(i,1)))/params.HBstops{i}{1,occIdx(i,1)};
end
% Proximity of primary hb day 1 versus primary hb day 2
row=1;
for i=1:2:size(params,1)
    primaryHBdist(row,1)=sqrt((params.HBcenter{i,1}{1,occIdx(i,1)}(1,1)-params.HBcenter{i+1,1}{1,occIdx(i+1,1)}(1,1))^2+...
        (params.HBcenter{i,1}{1,occIdx(i,1)}(1,2)-params.HBcenter{i+1,1}{1,occIdx(i+1,1)}(1,2))^2);
    row=row+1;
end

timeInCue=cell2mat(params.time_in_zone);

vars={'Subject','day','group','NumHBs','Occupancy','Entries',...
    'Distance2Cue','NumStops','AvgVelocity','area','avgStopDist',...
    'numCloseStops','clos2cue','time2hb','timeInCue','avgStopDuration'}; %varnames act as headers for hbMeas2Python
day=repmat({'D1';'D2'},size(params.HBcenter,1)/2,1); %create day variable 
group=[repmat({'Tg'},size(params.HBcenter,1)/2,1);...%create grouping variable
    repmat({'Wt'},size(params.HBcenter,1)/2,1)];
id=split(params.subID,'_'); %create subID as a within-subjects factor
[~,~,subject]=unique(id(:,1)); %create subID as a within-subjects factor 

hbMeas2Python=table(subject,day,group,numHB,primary_hbOcc,...
    primary_hbEntry, primary_hbDist2Cue,...
    primary_hbStops, primary_hbvel,primary_area,...
    primary_distHBstop,primary_numCloseHBstops,...
    primary_hbDistBinary,primary_time2HB,timeInCue,primary_hbStopDuration,'VariableNames',vars); %make table

writetable(hbMeas2Python,...
    'D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\hbData.csv'); %save data



%% Probe Trial analysis 

%difference scores

pathL = cell2mat(params.pathL); 
pathL_diff = pathL(2:2:end,:) - pathL(1:2:end,:);

searchArea = cell2mat(params.searchArea);
searchArea_diff = searchArea(2:2:end,:) - searchArea(1:2:end,:);

runSpeed = cell2mat(params.meanLinVel_run);
runSpeed_diff = runSpeed(2:2:end,:) - runSpeed(1:2:end,:);

runAcell = params.med_lin_accel; 
runAcell_diff = runAcell(2:2:end,:) - runAcell(1:2:end,:);

runAngVel = cell2mat(params.meanAngVel_run);
runAngVel_diff = runAngVel(2:2:end,:) - runAngVel(1:2:end,:);

stopAngVel = cell2mat(params.meanAngVel_stop);
stopAngVel_diff = stopAngVel(2:2:end,:) - stopAngVel(1:2:end,:);

numRuns = cell2mat(params.NumRuns);
numRuns_diff = numRuns(2:2:end,:) - numRuns(1:2:end,:);

numStops = cell2mat(params.NumStops); 
numStops_diff = numStops(2:2:end,:) - numStops(1:2:end,:);

outside_hb_stop = cell2mat(params.outside_hb_stops);
outside_hb_stop_diff = outside_hb_stop(2:2:end,:) - outside_hb_stop(1:2:end,:);

CueEntries = params.CueEntries;
CueEntries_diff = CueEntries(2:2:end,:) - CueEntries(1:2:end,:);

CueStops = params.CueStops;
CueStops_diff = CueStops(2:2:end,:) - CueStops(1:2:end,:);

InZoneTime= cell2mat(params.time_in_zone);
InZoneTime_diff = InZoneTime(2:2:end,:) - InZoneTime(1:2:end,:);

vars={'Subject','day','group','pathL_diff','searchArea_diff','runSpeed_diff','runAngVel_diff','stopAngVel_diff'...
    ,'numRuns_diff','numStops_diff','outside_hb_stop_diff','CueEntries_diff','CueStops_diff','runAcell_diff','cue_dwell_diff'};
%Now we'll create variables for levels of independent variables as well as
%make the unique subject id. 
day=repmat({'D2'},size(CueStops_diff,1),1); %create day variable 
group=[repmat({'Tg'},size(CueStops_diff,1)/2,1);...%create grouping variable
    repmat({'Wt'},size(CueStops_diff,1)/2,1)];
id=split(params.subID,'_'); %create subID as a within-subjects factor
subject=unique(id(:,1)); %create subID as a within-subjects factor 

%We can now create a table with only our new vars 
wholeTrial2Python=table(subject,day,group,pathL_diff,searchArea_diff,...
    runSpeed_diff,runAngVel_diff,stopAngVel_diff,numRuns_diff,numStops_diff,outside_hb_stop_diff,...
    CueEntries_diff,CueStops_diff,runAcell_diff,InZoneTime_diff,'VariableNames',vars); %make table

writetable(wholeTrial2Python,...
    'D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\wholeTrial_difference_day2.csv'); %save data


%Get logical of whether moving slow critera was met

for i=1:size(params.HBcenter,1)
    for ii = 1:size(params.hbOcc{i},2)
        move_slow{i,1}(1,ii) = cell2mat(params.HBclass{i}(1,ii)) > .75;
    end
end

%num of HB that meet above 
for i = 1 : length (params.hbOcc)
    numHB(i,1)= sum(move_slow{i});
end

for i=1:size(params.HBcenter,1)
    [~,occIdx(i,1)]=max(params.hbOcc{i}(move_slow{i})); %max given moving slow most of time;
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
    primary_hbStopDuration(i,1) = cell2mat(params.slowInHB{i}(1,occIdx(i,1)))/params.HBstops{i}{1,occIdx(i,1)};
end

timeInCue=cell2mat(params.time_in_zone);

vars={'Subject','day','group','NumHBs','Occupancy','Entries',...
    'Distance2Cue','NumStops','AvgVelocity','area','avgStopDist',...
    'numCloseStops','clos2cue','time2hb','timeInCue','avgStopDuration',}; %varnames act as headers for hbMeas2Python
day=repmat({'D1';'D2'},size(params.HBcenter,1)/2,1); %create day variable 
group=[repmat({'Tg'},size(params.HBcenter,1)/2,1);...%create grouping variable
    repmat({'Wt'},size(params.HBcenter,1)/2,1)];
id=split(params.subID,'_'); %create subID as a within-subjects factor
[~,~,subject]=unique(id(:,1)); %create subID as a within-subjects factor 

hbMeas2Python=table(subject,day,group,numHB,primary_hbOcc,...
    primary_hbEntry, primary_hbDist2Cue,...
    primary_hbStops, primary_hbvel,primary_area,...
    primary_distHBstop,primary_numCloseHBstops,...
    primary_hbDistBinary,primary_time2HB,timeInCue,primary_hbStopDuration,'VariableNames',vars); %make table

writetable(hbMeas2Python,...
    'D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\hbData.csv'); %save data

% SECOND PARAGRAPH RESULTS (PROBE TEST: PROXIMAL CUE REMOVAL). 
for i = 1:length(param_idx)
    temp = params.tsStop_cue{i};
    temp(isnan(params.tsStop_cue{i}),:)=[];
    tsStop_cue{i,1} = temp;
    clear temp
end
Tg_bin_entry = [];
Wt_bin_entry = [];
for i = 1:length(param_idx)
    temp_ts_entry(i,1) = params.tsEntry_cue{i}(1,1);
    temp_bin_entry(i,:) = params.bin_entries_cue{i};
end
% SECOND PARAGRAPH RESULTS (PROBE TEST: PROXIMAL CUE REMOVAL). 

day2_Tg = temp_ts_entry(contains(param_idx,{'Tg'}) & contains(param_idx,{'D2'}),:);
day2_Wt = temp_ts_entry(contains(param_idx,{'WT'}) & contains(param_idx,{'D2'}),:);

cue_measure = table(subject,day,group,temp_ts_entry,'VariableNames',{'subject','day','group','time2cue'});
writetable(cue_measure,...
    'D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\cueData.csv'); %save data

% SECOND PARAGRAPH RESULTS (PROBE TEST: PROXIMAL CUE REMOVAL). 
edges = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31];
edges = edges*60;
for i=1:length(param_idx)
   
    bin_dur_mat(i,:)=histcounts(params.tsEntry_cue{i},edges); 
    
end

prop_cue_day2_Tg = sum(sum(bin_dur_mat(contains(param_idx,{'Tg'}) & contains(param_idx,{'D2'}),:),2)~=0) ;
prop_cue_day2_Wt = sum(sum(bin_dur_mat(contains(param_idx,{'WT'}) & contains(param_idx,{'D2'}),:),2)~=0) ;

[h,p, chi2stat,df] = prop_test([prop_cue_day2_Tg prop_cue_day2_Wt],[12 12],'true')

% Proximity of primary hb day 1 versus primary hb day 2 (RESULTS LAST
% PARAGRAPH)
row=1;
for i=1:2:size(params,1)
    primaryHBdist(row,1)=sqrt((params.HBcenter{i,1}{1,occIdx(i,1)}(1,1)-params.HBcenter{i+1,1}{1,occIdx(i+1,1)}(1,1))^2+...
        (params.HBcenter{i,1}{1,occIdx(i,1)}(1,2)-params.HBcenter{i+1,1}{1,occIdx(i+1,1)}(1,2))^2);
    row=row+1;
end

subject = unique(subject,'sorted');
group=[repmat({'Tg'},size(primaryHBdist,1)/2,1);...%create grouping variable
    repmat({'Wt'},size(primaryHBdist,1)/2,1)];
cue_measure = table(subject,group,primaryHBdist,'VariableNames',{'subject','group','primaryHBdist'});
writetable(cue_measure,...
    'D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\HBdist.csv'); %save data
