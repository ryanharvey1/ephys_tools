
%% This script organizing data for the TgF344-AD open field study conducted by Berkowitz et al. 
%  Most data are exported as csv files for further analysis in Python.
%  However, some measures are directly analyzed below. 

load('D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\params_V18');
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
vars={'Subject','day','group','pathL','searchArea','runSpeed','numStops','CueEntries','CueStops'};

%Lets pull the variables from the table. 
pathL=cell2mat(params.pathL); 
searchArea=cell2mat(params.searchArea);
runSpeed=cell2mat(params.runSpeed);
numStops=cell2mat(params.NumStops); 
CueEntries=params.CueEntries;
CueStops=params.CueStops;

%Now we'll create variables for levels of independent variables as well as
%make the unique subject id. 
day=repmat({'D1';'D2'},size(params.time_in_zone,1)/2,1); %create day variable 
group=[repmat({'Tg'},size(params.time_in_zone,1)/2,1);...%create grouping variable
    repmat({'Wt'},size(params.time_in_zone,1)/2,1)];
id=split(params.subID,'_'); %create subID as a within-subjects factor
[~,~,subject]=unique(id(:,1)); %create subID as a within-subjects factor 

%We can now create a table with only our new vars 
wholeTrial2Python=table(subject,day,group,pathL,searchArea,runSpeed,numStops,CueEntries,CueStops,'VariableNames',vars); %make table

writetable(wholeTrial2Python,...
    'D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\wholeTrial_measures.csv'); %save data

%% Get Binned Stop time duration (#2)
edges = [0 2 10 60 1800];

for i=1:length(param_idx)
   bin_dur_mat(i,:)=histcounts(cell2mat(params.timeStopped{i})',edges); 
    
end

% save vars as table 

vars_continuous = {'Subject','day','group','dur_bin','stop_duration'};

dur_bin=[zeros(48,1)+1; zeros(48,1)+2; zeros(48,1)+3; zeros(48,1)+4];

continuous_bin_dur = table([subject;subject;subject;subject], [day;day;day;day],...
    [group;group;group;group], dur_bin, ...
    [bin_dur_mat(:,1); bin_dur_mat(:,2); bin_dur_mat(:,3); bin_dur_mat(:,4)]...
    ,'VariableNames',vars_continuous);

writetable(continuous_bin_dur,...
    'D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\continuous_bin_dur.csv'); %save data

varargout=shadedErrorBar(x,y,errBar,lineProps,transparent)

fig = figure; 
fig.Color = [1 1 1];
X = categorical({'0-2','2-10','10-60','60+'});
X = reordercats(X,{'0-2','2-10','10-60','60+'});

subplot(2,1,1)
b=bar(X,[norm_wt1_bin' norm_tg1_bin'])
b(2).FaceColor =  'r';
b(1).FaceColor = [.25 .25 .25];
box off
xlabel('Time (seconds)')
ylabel('Proportion of Stops')
title('Day 1')
legend({'WT','TgF344-AD'}); legend('boxoff')
set(gca,'FontSize',14,'FontName','Helvetica','FontWeight','bold','LineWidth',2)
subplot(2,1,2)
b=bar(X,[norm_wt2_bin' norm_tg2_bin'])
b(2).FaceColor =  'r';
b(1).FaceColor = [.25 .25 .25];
box off
xlabel('Time (seconds)')
ylabel('Proportion of Stops')
title('Day 2')
legend({'WT','TgF344-AD'});legend('boxoff')
set(gca,'FontSize',14,'FontName','Helvetica','FontWeight','bold','LineWidth',2)

 export_fig('D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Figures\PNGs\Binned_Stops.png','-m4') 
 
% timeStopped_tg1=horzcat(params.timeStopped{contains(param_idx,'Tg') & contains(param_idx,'D1')});
% timeStopped_wt1=horzcat(params.timeStopped{contains(param_idx,'WT') & contains(param_idx,'D1')});
% timeStopped_tg2=horzcat(params.timeStopped{contains(param_idx,'Tg') & contains(param_idx,'D2')});
% timeStopped_wt2=horzcat(params.timeStopped{contains(param_idx,'WT') & contains(param_idx,'D2')});
% 
% tg1_bin = histcounts(cell2mat(timeStopped_tg1)',edges);
% norm_tg1_bin = tg1_bin/size(cell2mat(timeStopped_tg1)',1);
% 
% wt1_bin = histcounts(cell2mat(timeStopped_wt1)',edges);
% norm_wt1_bin = wt1_bin/size(cell2mat(timeStopped_wt1)',1);
% 
% tg2_bin = histcounts(cell2mat(timeStopped_tg2)',edges);
% norm_tg2_bin = tg2_bin/size(cell2mat(timeStopped_tg2)',1);
% 
% wt2_bin = histcounts(cell2mat(timeStopped_wt2)',edges);
% norm_wt2_bin = wt2_bin/size(cell2mat(timeStopped_wt2)',1);
%  data_cue = [tg1_bin; wt1_bin];
% 
% contingency_cue = table(subject, day, group, data_cue(:,1),data_cue(:,2),data_cue(:,3),data_cue(:,4)...
%     ,'VariableNames',vars,'RowNames',{'TgF344-AD','F344'});
% 
% data_probe = [tg2_bin; wt2_bin];
% 
% contingency_probe = table(subject, day, group,  data_probe(:,1),data_probe(:,2),data_probe(:,3),data_probe(:,4)...
%     ,'VariableNames',vars,'RowNames',{'TgF344-AD','F344'});
% writetable(contingency_cue,...
%     'D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\binned_duration_cue.csv'); %save data
% writetable(contingency_probe,...
%     'D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\binned_duration_probe.csv'); %save data


%% Get Primary home base measures lb 9/18/19 Add mean stop duration in home base (#3) 
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
    primary_hbStopDuration(i,1) = cell2mat(params.slowInHB{i}(1,occIdx(i,1)))/params.HBstops{i}{1,occIdx(i,1)};
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
    'numCloseStops','clos2cue','time2hb','timeInCue','avtStopDuration'}; %varnames act as headers for hbMeas2Python
day=repmat({'D1';'D2'},size(params.HBcenter,1)/2,1); %create day variable 
group=[repmat({'Tg'},size(params.HBcenter,1)/2,1);...%create grouping variable
    repmat({'Wt'},size(params.HBcenter,1)/2,1)];
id=split(params.subID,'_'); %create subID as a within-subjects factor
[~,~,subject]=unique(id(:,1)); %create subID as a within-subjects factor 

hbMeas2Python=table(subject,day,group,numHb,primary_hbOcc,...
    primary_hbEntry, primary_hbDist2Cue,...
    primary_hbStops, primary_hbvel,primary_area,...
    primary_distHBstop,primary_numCloseHBstops,...
    primary_hbDistBinary,primary_time2HB,timeInCue,primary_hbStopDuration,'VariableNames',vars); %make table

writetable(hbMeas2Python,...
    'D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\hbData.csv'); %save data

%% get velocity measures during stopping (item 5)


%standard error of linear and angular velocity
param_idx=params.subID;
ang_vel_Tg=[];
ang_vel_WT=[];
for i=1:length(param_idx)
    runIdx=params.runIdx{i};
    temp_angVel=params.angVel{i};
    temp_linVel=params.pathIV{i};
    
    if contains(param_idx(i),'Tg')
        motion(i,1)=1;
        ang_vel_Tg=[ang_vel_Tg; temp_angVel];
    else
        motion(i,1)=0;
        ang_vel_WT=[ang_vel_WT; temp_angVel];
    end
    
    if mod(i,2)==0
        motion(i,2)=2;
    else
        motion(i,2)=1;
    end
    
    motion(i,3)=nanstd(abs(temp_angVel));
    motion(i,4)=nanstd(abs(temp_linVel));
    motion(i,5)=nanstd(abs(temp_angVel(runIdx)));
    motion(i,6)=nanstd(abs(temp_linVel(runIdx)));
    motion(i,7)=nanmean(abs(temp_angVel(~runIdx)));
    motion(i,8)=nanstd(abs(temp_angVel(~runIdx)));
    
end

vars={'group','day','mean_abs_vel','std_abs_vel'};

velocity_data = table(motion(:,1),motion(:,2),...
    motion(:,7),motion(:,8),'VariableNames',vars);

writetable(velocity_data,...
    'D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\velocity_data.csv'); %save data


%% Number of Rats with x num of home bases 

% Get home bases that are minimum 5 min in duration 

for i = 1 : length (params.hbOcc)

    numHB(i,1)= sum(params.hbOcc{i} > 300);
    
end

day_idx = zeros(48,1);
day_idx([2:2:48]',1)=1; 
group_idx = [ones(24,1) ;zeros(24,1)];

vars = {'group','day','number_of_homebase'};
num_hb = table(group_idx, day_idx, numHB,'VariableNames',vars);

writetable(num_hb,...
    'D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\num_hb.csv'); %save data


edges = [1 2 3 4];
tg_binned_1 = histcounts(numHB(group_idx == 1 & day_idx == 0,1),edges);
wt_binned_1 = histcounts(numHB(group_idx == 0 & day_idx == 0,1),edges);
tg_binned_2 = histcounts(numHB(group_idx == 1 & day_idx == 1,1),edges);
wt_binned_2 = histcounts(numHB(group_idx == 0 & day_idx == 1,1),edges);




fig = figure; 
fig.Color = [1 1 1];

labels = {'One','Two','Three'};
subplot(2,2,1)
tg_d1 = numHB(day_idx == 0 & group_idx ==1,1);
pie(histcounts(tg_d1,[1 2 3 4]))
legend(labels,'Orientation','horizontal','Location','southoutside')
legend('boxoff')
title('Tg-F344-AD Day 1')
set(gca,'FontSize',13,'FontName','Helvetica','FontWeight','bold')
subplot(2,2,2)
wt_d1 = numHB(day_idx == 0 & group_idx ==0,1);
pie(histcounts(wt_d1,[1 2 3 4]))
legend(labels,'Orientation','horizontal','Location','southoutside')
legend('boxoff')
title('F344 Day 1')
set(gca,'FontSize',13,'FontName','Helvetica','FontWeight','bold')
subplot(2,2,3)
tg_d2 = numHB(day_idx == 1 & group_idx ==1,1);
pie(histcounts(tg_d2,[1 2 3 4]))
legend(labels,'Orientation','horizontal','Location','southoutside')
legend('boxoff')
title('Tg-F344-AD Day 2')
set(gca,'FontSize',13,'FontName','Helvetica','FontWeight','bold')
subplot(2,2,4)
wt_d2 = numHB(day_idx == 1 & group_idx ==0,1);
pie(histcounts(wt_d2,[1 2 3 4]))
legend(labels,'Orientation','horizontal','Location','southoutside')
legend('boxoff')
title('F344 Day 2')
set(gca,'FontSize',13,'FontName','Helvetica','FontWeight','bold')

export_fig('D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Figures\PNGs\Home_base_pie.png','-m4') 

%% Probe Trial analysis 
CueEntries_tg = vertcat(params.CueEntries(contains(param_idx,'Tg') & contains(param_idx,'D2')));
CueEntries_wt = vertcat(params.CueEntries(contains(param_idx,'WT') & contains(param_idx,'D2')));

CueStops_tg = vertcat(params.CueStops(contains(param_idx,'Tg') & contains(param_idx,'D2')));
CueStops_wt = vertcat(params.CueStops(contains(param_idx,'WT') & contains(param_idx,'D2')));

BinStops = vertcat(params.bin_stops{contains(param_idx,'D2')});

BinEntries = vertcat(params.bin_entries{contains(param_idx,'D2')});

vars={'Subject','group','stops'};

subject=[1:24]';
group = [ones(12,1);zeros(12,1)];

bin_stop=table(subject,group,BinStops,'VariableNames',vars);

vars={'Subject','group','entries'};

bin_entries=table(subject,group,BinEntries(:,1),'VariableNames',vars);

writetable(bin_stop,...
    'D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\bin_stop.csv'); %save data

writetable(bin_entries,...
    'D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\bin_entry.csv'); %save data



fig = figure; 
fig.Color=[1 1 1];

subplot(2,2,1)
plotspread_wrapper(CueEntries_wt,CueEntries_tg,{'F344','TgF344-AD'})
set(gca,'FontWeight','bold','FontName','Helvetica','FontSize',14)
ylim([0 max([CueEntries_wt;CueEntries_tg;CueStops_wt;CueStops_tg])])
title('Entries')

subplot(2,2,2)
plotspread_wrapper(CueStops_wt,CueStops_tg,{'F344','TgF344-AD'})
set(gca,'FontWeight','bold','FontName','Helvetica','FontSize',14)
ylim([0 max([CueEntries_wt;CueEntries_tg;CueStops_wt;CueStops_tg])])
title('Stops')

subplot(2,2,3)
plotspread_wrapper(BinEntries(12:end,:),BinEntries(1:12,:),{'F344','TgF344-AD'})
set(gca,'FontWeight','bold','FontName','Helvetica','FontSize',14)
ylim([0 max([BinEntries;BinStops])])
title('Entries first 2 minutes')

subplot(2,2,4)
plotspread_wrapper(BinStops(12:end,:),BinStops(1:12,:),{'F344','TgF344-AD'})
set(gca,'FontWeight','bold','FontName','Helvetica','FontSize',14)
ylim([0 max([BinEntries;BinStops])])
title('Stops first 2 minutes')

export_fig('D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Figures\PNGs\entries_stops.png','-m4') 


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

% figure; 
% subplot(1,2,1)
% plotspread_wrapper(cell2mat(wt1),cell2mat(tg1),{'WT','Tg'})
% title('Time to first stop in Primary HB - Cue Present')
% ylim([0 max([cell2mat(wt1);cell2mat(tg1);cell2mat(wt2);cell2mat(tg2)])])
% ylabel('time (s)')
% subplot(1,2,2)
% plotspread_wrapper(cell2mat(wt2),cell2mat(tg2),{'WT','Tg'})
% title('Time to first stop in Primary HB - Cue Absent')
% ylim([0 max([cell2mat(wt1);cell2mat(tg1);cell2mat(wt2);cell2mat(tg2)])])
% ylabel('time (s)')
% 
% 
fig = figure; 
fig.Color = [1 1 1];
subplot(1,2,1)
plotspread_wrapper(primary_hbStopDuration(25:2:end,1),primary_hbStopDuration(1:2:24,1),{'WT','Tg'})
title('Day 1')
ylabel('Average Stop Duration in Primary Home Base (s)')
ylim([0 (max(primary_hbStopDuration)+10)])
set(gca,'LineWidth',2,'FontSize',14,'FontName','Helvetica')
subplot(1,2,2)
plotspread_wrapper(primary_hbStopDuration(26:2:end,1),primary_hbStopDuration(2:2:24,1),{'WT','Tg'})
title('Day 2')
legend({'WT','TgF344-AD'});
legend('boxoff')
ylim([0 (max(primary_hbStopDuration)+10)])
set(gca,'LineWidth',2,'FontSize',14,'FontName','Helvetica')
export_fig('D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Figures\PNGs\AveStop_inHB.png','-m4') 

fig=figure; 
fig.Color=[1 1 1];
scatter(primary_hbStopDuration(1:2:24,1),primary_hbStopDuration(2:2:24,1),'r','filled')
hold on; 
scatter(primary_hbStopDuration(25:2:end,1),primary_hbStopDuration(26:2:end,1),'k','filled')


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
% 
% subplot(1,2,1)
% plotspread_wrapper(sum(wt1(:,5:7),2),sum(tg1(:,5:7),2),{'WT','Tg'})
% title('Cue quadrant dwell time')
% 
% subplot(1,2,2)
% plotspread_wrapper(sum(wt2(:,5:7),2),sum(tg2(:,5:7),2),{'WT','Tg'})
% title('Cue quadrant dwell time')

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
